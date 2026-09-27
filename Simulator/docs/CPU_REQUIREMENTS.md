# CPU requirements

The standard C++ build targets `x86-64-v2`, which runs on any Intel or AMD x86-64 machine from about 2009 onwards. The engine needs that level for one instruction only: `POPCNT`, which counts the set bits in a word and is used on every event to count the allowed bonds.  The code still runs on any non-x86 machine and gives the same answers there -- the compiler simply replaces each hardware count with a software loop over the bits, so only the speed becomes much slower. See the end of this file for how.

## BMI2, chosen at run time rather than at build time

Picking the `n`-th allowed bond out of a mask word is the step that dominates the cost of an event. There are two ways to do it.

**The portable version** clears the lowest set bit `n` times, then reports the position of whichever bit is now lowest:

```cpp
while (n-- > 0) mask &= mask - 1;      // mask & (mask-1) clears the lowest set bit
return std::countr_zero(mask);         // position of the lowest remaining set bit
```

It is ordinary C++ and compiles anywhere, but each iteration depends on the result of the previous one, so the cost grows with `n`: finding the `n`-set bit takes `n` steps.

**The `PDEP` version** does the whole search in one machine instruction. `PDEP` ("parallel bit deposit") takes two words. It walks the *mask* from the low end, and every time it meets a set bit, it consumes the next bit of the *source* and writes it in that position; where the mask is zero the output is zero. So the source's bits are spread out into the mask's set positions, in order:

```text
mask     1001011010000000     set at positions 7, 9, 10, 12, 15
source   0000000000001011     its low bits, from the bottom: 1, 1, 0, 1
PDEP     0001001010000000     those bits land at 7, 9, 10, 12
```

The selection is the special case where the source holds a single `1` at position `n`. There is then only one bit to deposit, and it lands in the `n`-th set position of the mask, so the result is a word with exactly one bit set, at the position we were looking for:

```text
n = 0    1<<0 = ...0001   ->  0000000010000000   position  7   the 1st allowed bond
n = 1    1<<1 = ...0010   ->  0000001000000000   position  9   the 2nd
n = 2    1<<2 = ...0100   ->  0000010000000000   position 10   the 3rd
n = 3    1<<3 = ...1000   ->  0001000000000000   position 12   the 4th
```

`std::countr_zero` then turns that one word into the index itself by counting the zeros below the surviving bit:

```cpp
return std::countr_zero(_pdep_u64(std::uint64_t{1} << n, mask));
```

The hardware performs the scatter for all 64 positions at once — hence "parallel" — rather than visiting the set bits one by one, so its cost does not depend on `n` at all.

## Why `PDEP` needs the attribute

`PDEP` belongs to the BMI2 instruction set, which arrives only at `x86-64-v3`. The `x86-64-v2` baseline this project builds against provides `POPCNT` but nothing beyond it. So the compiler refuses to emit `PDEP` in ordinary code (even when the CPU compiling it has BMI2), and a build that simply used it would fail with `target specific option mismatch`.

`__attribute__((target("bmi2")))` lifts that restriction for the single function it applies to — the four-line `select_in_word_pdep` — while everything else stays at the `x86-64-v2` baseline. It grants permission to *compile* the instruction and says nothing about the machine that will later run it — which is harmless, because an instruction only matters once it is executed. On a CPU without BMI2, there is a run-time check that prevents the code from ever taking that branch, the `PDEP` bytes are never fetched, and the portable walk runs instead: the binary is still perfectly usable, it simply takes the slower route. 

## Why the choice is made at run time

Two questions have to be answered before `PDEP` can be used, and only the first can be settled while compiling. **Does this build contain the instruction?** — yes, thanks to the attribute. **Should this machine execute it?** — that depends on the processor in front of you, which the compiler cannot know.

The second question is not simply "does the CPU support BMI2", because support does not imply speed. `PDEP` is a single fast instruction on Intel from Haswell (2013) and on AMD from Zen 3 (2020), but on Zen 1 and Zen 2 it is microcoded and takes hundreds of cycles -- far slower than the loop it would replace. So a CPU can advertise BMI2 and still be better off with the portable walk.

`detect_usable_pdep()` therefore asks for BMI2 support, then reads the vendor and family from `CPUID` and rejects AMD family `0x17`, which is Zen 1 and Zen 2. It runs once, when the extension module is loaded, and stores its verdict in a boolean. Every event then reads that boolean and branches.

## Measured

Selecting one bit from a 64-bit word, averaged over 13 million calls on an Intel i7-8550U. Both paths select the same bond for the same draw, so a run is bitwise identical either way; only the speed differs. The two are compared at three fillings and three choices of `n`, since that is what decides the outcome:

| bits set in the word | which set bit  | portable  | `PDEP`   | speed-up  |
| -------------------- | -------------- | --------- | -------- | --------- |
| 8 of 64              | the 1st        | 8.97 ns   | 12.28 ns | 0.7 times |
| 8 of 64              | the middle one | 12.84 ns  | 15.98 ns | 0.8 times |
| 8 of 64              | the last       | 27.22 ns  | 13.90 ns | 2.0 times |
| 32 of 64             | the 1st        | 6.08 ns   | 12.50 ns | 0.5 times |
| 32 of 64             | the middle one | 42.23 ns  | 20.07 ns | 2.1 times |
| 32 of 64             | the last       | 105.24 ns | 23.31 ns | 4.5 times |
| 56 of 64             | the 1st        | 10.49 ns  | 20.96 ns | 0.5 times |
| 56 of 64             | the middle one | 79.11 ns  | 16.09 ns | 4.9 times |
| 56 of 64             | the last       | 158.75 ns | 17.91 ns | 8.9 times |

The pattern is the one the two algorithms predict: `PDEP` is flat, the portable walk is not. For the first set bit the loop wins, because it exits immediately and `PDEP` still has to issue. From about the middle onwards `PDEP` is ahead, and the denser the word the further ahead it gets. In a typical run roughly half the sites are occupied and the bond is drawn uniformly, so the one being looked for usually sits somewhere in the middle of the word -- the region where `PDEP` is ahead. Over a whole simulation that is worth about 15-30% in speed.

## Pre-historic computer (non-x86 computer)

An old x86-64 CPU without x86-64-v2 support will abort with an illegal instruction. A non-x86 computer, such as an ARM laptop, may reject the compiler flag. In either case, open `setup_cpp.py` and change

```python
extra_compile_args=['-O3', '-std=c++20', '-march=x86-64-v2'],
```

to

```python
extra_compile_args=['-O3', '-std=c++20'],
```

Then rebuild normally. No C++ changes are required: `std::popcount` uses the best implementation available for the selected compiler target. 
