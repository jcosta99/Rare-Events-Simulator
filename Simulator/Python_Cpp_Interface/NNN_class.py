"""Thin workflow wrapper around the complete C++ scientific engine."""

import base64
import csv
import json
from pathlib import Path

try:
    from . import _engine
except ImportError as error:
    raise ImportError(
        'The C++ engine is not built. Run: '
        'conda run -n scientific-python python setup_cpp.py build_ext --inplace'
    ) from error


class NNN:
    """Run the rejection-free C++ implementation through a Python interface.

    The output-row schema is documented in ``docs/OUTPUT_CONVENTIONS.md``, the
    physical parameters in ``theory_algorithm_implementation.tex``.
    """

    def __init__(self, L=20, M=100, tmax=20.0, E=0.0, s=0.0, k=0.0,
                 model=None, epsilon=None, delta_kls=None,
                 right_weights=None, left_weights=None,
                 filling=0.5,
                 alpha=1.0, gamma=1.0, delta=1.0, beta=1.0,
                 target_min=None, target_max=None,
                 cloning_interval=None, record_interval=None,
                 seed=12345, initial='alternating', progress=None,
                 dynamics='direct',
                 boundary='open'):
        if L < 2:
            raise ValueError('L must be at least 2')
        if M < 2:
            raise ValueError('M must be at least 2')
        if tmax <= 0:
            raise ValueError('tmax must be positive')
        if not 0 <= filling <= 1:
            raise ValueError('filling must lie between 0 and 1')
        if min(k, alpha, gamma, delta, beta) < 0:
            raise ValueError('k and boundary rates must be non-negative')
        # Explicit weights, when given, replace the KLS formula entirely,
        # so the epsilon/delta bounds do not apply to them.
        explicit = right_weights is not None or left_weights is not None
        # The model names the rate family and, with it, the bookkeeping the
        # engine uses.  Leaving it unset keeps the older behaviour of reading
        # the family off the rates themselves; 'WASEP' has to be asked for,
        # because it is a different code path and not just epsilon = 0.
        if model is None:
            model = 'NNN' if explicit else 'KLS'
        model = str(model).upper()
        if model not in {'WASEP', 'KLS', 'NNN'}:
            raise ValueError("model must be 'WASEP', 'KLS' or 'NNN'")
        # Unset interaction parameters mean the historical KLS defaults, but
        # only for a KLS run: a WASEP or NNN run must not silently acquire an
        # interaction it does not use.
        if model == 'KLS':
            if epsilon is None:
                epsilon = 0.6
            if delta_kls is None:
                delta_kls = 0.0
        if model == 'WASEP':
            if explicit:
                raise ValueError('a WASEP run has one rate per direction; '
                                 'do not pass right_weights or left_weights')
            if epsilon or delta_kls:
                raise ValueError('a WASEP run has no interaction; use '
                                 "model='KLS' for epsilon or delta_kls")
            # The eight class weights are all one, so nothing else is needed.
            epsilon = delta_kls = 0.0
        elif model == 'NNN':
            if not explicit:
                raise ValueError(
                    "model='NNN' needs right_weights and left_weights")
            if epsilon or delta_kls:
                raise ValueError('explicit weights replace the KLS form; '
                                 'do not pass epsilon or delta_kls with them')
            epsilon = delta_kls = 0.0
        elif explicit:
            raise ValueError("model='KLS' is set by epsilon and delta_kls; "
                             "pass model='NNN' to give weights directly")
        if explicit:
            if right_weights is None or left_weights is None:
                raise ValueError('give right_weights and left_weights together')
            right_weights = [float(w) for w in right_weights]
            left_weights = [float(w) for w in left_weights]
            if len(right_weights) != 4 or len(left_weights) != 4:
                raise ValueError('right_weights and left_weights need four '
                                 'rates each, ordered (0,0) (0,1) (1,0) (1,1)')
            if min(right_weights + left_weights) <= 0:
                raise ValueError('every hop rate must be positive')
        # The four KLS weights are 1 + delta, 1 - epsilon, 1 + epsilon and
        # 1 - delta: a nonzero (a - d) forces a + d = 1, which kills the delta
        # term, so the two parameters never appear in the same weight.  Each
        # therefore only has to stay within (-1, 1) on its own; requiring
        # |epsilon| + |delta| < 1 as well would reject valid models such as
        # epsilon = delta = 0.6, whose rates are 1.6, 0.4, 1.6, 0.4.
        elif not -1.0 < epsilon < 1.0:
            raise ValueError('epsilon must lie strictly between -1 and 1 '
                             'so that every hop rate is positive')
        elif not -1.0 < delta_kls < 1.0:
            raise ValueError('delta_kls must lie strictly between -1 and 1 '
                             'so that every hop rate is positive')
        # Two weight-ratio targets: resample once the population's weights
        # differ by target_min, never let them differ by more than target_max.
        # Supplying only target_max implies target_min as its two-thirds power.
        if target_max is None:
            target_max = 10.0
        if target_min is None:
            target_min = target_max ** (2.0 / 3.0)
        if target_min <= 1:
            raise ValueError('target_min must be greater than 1')
        if target_max <= target_min:
            raise ValueError('target_max must be greater than target_min')
        if cloning_interval is not None and cloning_interval <= 0:
            raise ValueError('cloning_interval must be positive')
        if record_interval is not None and record_interval <= 0:
            raise ValueError('record_interval must be positive')
        if progress is not None and (
                isinstance(progress, bool)
                or not isinstance(progress, int)
                or progress < 1):
            raise ValueError('progress must be a positive integer or None')
        if initial not in {'empty', 'alternating', 'random', 'full'}:
            raise ValueError(f'Unknown initial condition: {initial}')
        if dynamics not in {'heap', 'uniformized', 'direct'}:
            raise ValueError('dynamics must be heap, uniformized or direct')
        if boundary not in {'open', 'periodic'}:
            raise ValueError('boundary must be open or periodic')
        self.parameters = {
            'L': L, 'M': M, 'tmax': tmax, 'E': E, 's': s, 'k': k,
            # 'WASEP', 'KLS' or 'NNN'; see the constructor for what each means.
            'model': model,
            # KLS interaction; both zero reduces the model to WASEP.
            'epsilon': epsilon, 'delta_kls': delta_kls,
            'right_weights': right_weights, 'left_weights': left_weights,
            'filling': filling,
            'alpha': alpha, 'gamma': gamma, 'delta': delta, 'beta': beta,
            # Resample once the heaviest walker weight exceeds the lightest by
            # target_min; never let that ratio pass target_max.
            'target_min': target_min,
            'target_max': target_max,
            # Accepted for API compatibility. Neither direct C++ event loop
            # needs NumPy-style bounded proposal batches.
            'cloning_interval': cloning_interval,
            'record_interval': record_interval,
            'seed': seed, 'initial': initial, 'progress': progress,
            'dynamics': dynamics,
            'boundary': boundary,
        }
        self._simulation = _engine.create(self.parameters)
        self.saved_rows = None
        self.cache_action = 'new simulation'
        self.source_csv_path = None
        self.source_json_path = None

    def run(self):
        """Run or resume to ``tmax`` and return the complete time series."""

        new_rows = _engine.run(self._simulation)
        for row in new_rows:
            row['DensityProfile'] = json.dumps(row['DensityProfile'])
        if self.saved_rows is None:
            self.saved_rows = new_rows
        else:
            self.saved_rows.extend(new_rows)
        return list(self.saved_rows)

    def checkpoint(self):
        """Return opaque bytes containing the exact C++ restart state."""

        return _engine.checkpoint(self._simulation)

    def restore_checkpoint(self, checkpoint):
        """Restore an exact checkpoint produced by this C++ engine."""

        _engine.restore(self._simulation, checkpoint)

    def cache_parameters(self):
        """Return run-defining values except final time and progress output."""

        parameters = {
            name: value for name, value in self.parameters.items()
            if name not in {'tmax', 'progress'}
        }
        if self.parameters['boundary'] == 'periodic':
            # Reservoirs and the legacy initial selector have no dynamics on
            # a fixed-filling closed ring.
            for name in ('alpha', 'gamma', 'delta', 'beta', 'initial'):
                parameters.pop(name)
        else:
            # Filling only defines the initial particle number of a ring.
            parameters.pop('filling')
        # Only the rate description that is actually in force may take part in
        # the cache signature.  Explicit weights make epsilon and delta inert,
        # the KLS parameters leave the weight lists empty, and a WASEP run uses
        # neither, so keeping the unused ones would block reuse between
        # otherwise identical runs.
        model = self.parameters['model']
        if model != 'NNN':
            parameters.pop('right_weights')
            parameters.pop('left_weights')
        if model != 'KLS':
            parameters.pop('epsilon')
            parameters.pop('delta_kls')
        return parameters

    def constructor_parameters(self, tmax=None):
        """Return constructor arguments reproducing this engine."""

        parameters = dict(self.parameters)
        if tmax is not None:
            parameters['tmax'] = tmax
        return parameters

    def rate_family(self):
        """Return the model name: 'WASEP', 'KLS' or 'NNN'."""

        return self.parameters['model']

    def result_stem(self):
        """Build the canonical boundary-labelled NNN filename."""

        p = self.parameters
        boundary_code = 'C' if p['boundary'] == 'periodic' else 'O'
        # The model name is the filename family: WASEP for one rate per
        # direction, KLS for the Katz-Lebowitz-Spohn form, NNN for arbitrary
        # class weights.
        prefix = f"{self.rate_family()}_{boundary_code}"
        # The event mechanism is named too.  Only a run of the same mechanism
        # may be extended from a stored one -- their restart states are not
        # interchangeable -- so the name has to make the distinction visible
        # rather than leaving it buried in the metadata.
        prefix += f"_{p['dynamics'].capitalize()}"
        # A ring ignores ``initial`` altogether -- its particle number is fixed
        # by the filling -- so the tag belongs only to open chains, exactly as
        # in cache_parameters(). Naming it on a ring would give two identical
        # runs two different files.
        boundary_parameters = (
            f"_rho{p['filling']}" if p['boundary'] == 'periodic' else
            f"_alpha{p['alpha']}_gamma{p['gamma']}"
            f"_delta{p['delta']}_beta{p['beta']}_init{p['initial']}"
        )
        # The rate description goes into the filename so that runs of
        # different models never collide, and so that a cached result can be
        # recognised from its name alone.
        if p['model'] == 'NNN':
            right = '-'.join(f'{w:g}' for w in p['right_weights'])
            left = '-'.join(f'{w:g}' for w in p['left_weights'])
            rates = f'_R{right}_W{left}'
        elif p['model'] == 'KLS':
            rates = f"_eps{p['epsilon']}_del{p['delta_kls']}"
        else:
            # A WASEP run has no rates to name, which also keeps its filenames
            rates = ''
        return (
            f"{prefix}{rates}"
            f"_L{p['L']}_M{p['M']}_T{p['tmax']}"
            f"_E{p['E']}_s{p['s']}_k{p['k']}"
            f"{boundary_parameters}_seed{p['seed']}"
        )

    def result_paths(self, data_dir):
        # The stem ends in values such as "_rho0.5", whose trailing ".5"
        # Path.with_suffix would treat as an extension and replace, silently
        # truncating the density out of the filename.
        stem = self.result_stem()
        directory = Path(data_dir)
        return directory / f'{stem}.csv', directory / f'{stem}.json'

    def metadata(self):
        """Return flat, analysis-ready metadata and a restart checkpoint."""

        p = self.parameters
        diagnostics = self.diagnostics()
        closed = p['boundary'] == 'periodic'
        particles = int(p['filling'] * p['L'] + 0.5) if closed else None
        return {
            'algorithm': {
                'heap': 'C++ rejection-free accepted-event heap with '
                        'systematic resampling',
                'uniformized': 'C++ uniformized rejected proposals with '
                               'systematic resampling',
                'direct': 'C++ rejection-free per-walker evolution between '
                          'fixed cloning instants, with systematic resampling',
            }[p['dynamics']],
            'metadata_layout': 'flat-wasep-v1',
            'engine': f"nnn-{p['dynamics']}-v6",
            'boundary_condition': p['boundary'],
            'length': p['L'],
            'walkers': p['M'],
            'tmax': p['tmax'],
            'field_E': p['E'],
            'bias_s': p['s'],
            'measurement_strength_k': p['k'],
            'model': p['model'],
            'rate_model': diagnostics['rate_model'],
            'kls_epsilon': p['epsilon'] if p['model'] == 'KLS' else None,
            'kls_delta': p['delta_kls'] if p['model'] == 'KLS' else None,
            'right_weights': diagnostics['right_weights'],
            'left_weights': diagnostics['left_weights'],
            'filling': p['filling'] if closed else None,
            'particle_number': particles,
            'realized_filling': particles / p['L'] if closed else None,
            'alpha_left_injection': None if closed else p['alpha'],
            'gamma_left_removal': None if closed else p['gamma'],
            'delta_right_injection': None if closed else p['delta'],
            'beta_right_removal': None if closed else p['beta'],
            'target_min': p['target_min'],
            'target_max': p['target_max'],
            'record_interval': diagnostics['record_interval'],
            'seed': p['seed'],
            'initial_condition': (
                'uniform-fixed-particle-number' if closed else p['initial']
            ),
            'dynamics': p['dynamics'],
            'cache_parameters': self.cache_parameters(),
            'diagnostics': diagnostics,
            'cache_action': self.cache_action,
            'checkpoint': {
                'version': 10,
                'encoding': 'base64 opaque C++ binary',
                'data': base64.b64encode(self.checkpoint()).decode('ascii'),
            },
        }

    def write_results(self, csv_path, overwrite=False):
        """Write the current rows and restart state to CSV and JSON."""

        if self.saved_rows is None:
            raise RuntimeError('run() must be called before write_results()')
        csv_path = Path(csv_path)
        boundary_code = 'C' if self.parameters['boundary'] == 'periodic' else 'O'
        expected = f'{self.rate_family()}_{boundary_code}_'
        if not csv_path.name.startswith(expected):
            raise ValueError(f'Output filename must start with {expected}')
        json_path = csv_path.with_suffix('.json')
        if not overwrite and (csv_path.exists() or json_path.exists()):
            raise FileExistsError(f'Output already exists: {csv_path}')
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with csv_path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.saved_rows[0]))
            writer.writeheader()
            writer.writerows(self.saved_rows)
        with json_path.open('w') as stream:
            json.dump(self.metadata(), stream, indent=2, sort_keys=True)
            stream.write('\n')
        return csv_path, json_path

    @staticmethod
    def read_csv_rows(path):
        with Path(path).open(newline='') as stream:
            return list(csv.DictReader(stream))

    @classmethod
    def from_data_folder(cls, data_dir='data', **parameters):
        """Reuse the longest compatible C++ checkpoint, extending only time."""

        requested = cls(**parameters)
        signature = requested.cache_parameters()
        candidates = []
        data_dir = Path(data_dir)
        if data_dir.exists():
            # Results are named after the model: WASEP_[OC]_*, KLS_[OC]_* when
            # the rates come from the epsilon/delta form, and NNN_[OC]_* when
            # they were given explicitly. All three families are scanned; the
            # cache signature below is what actually decides compatibility, and
            # it includes both the model and the engine version.
            for json_path in sorted(set(data_dir.rglob('WASEP_*.json'))
                                    | set(data_dir.rglob('KLS_*.json'))
                                    | set(data_dir.rglob('NNN_*.json'))):
                try:
                    with json_path.open() as stream:
                        metadata = json.load(stream)
                except (OSError, ValueError, TypeError):
                    continue
                csv_path = json_path.with_suffix('.csv')
                if (
                    metadata.get('engine')
                    == f"nnn-{requested.parameters['dynamics']}-v6"
                    and metadata.get('cache_parameters') == signature
                    and metadata.get('checkpoint', {}).get('version') == 10
                    and csv_path.exists()
                ):
                    candidates.append((metadata, csv_path, json_path))
        if not candidates:
            return requested

        target_tmax = max(
            requested.parameters['tmax'],
            max(float(item[0]['tmax']) for item in candidates),
        )
        metadata, csv_path, json_path = max(
            (item for item in candidates if float(item[0]['tmax']) <= target_tmax),
            key=lambda item: float(item[0]['tmax']),
        )
        simulation = cls(**requested.constructor_parameters(tmax=target_tmax))
        checkpoint = base64.b64decode(metadata['checkpoint']['data'])
        simulation.restore_checkpoint(checkpoint)
        simulation.saved_rows = cls.read_csv_rows(csv_path)
        simulation.source_csv_path = csv_path
        simulation.source_json_path = json_path
        simulation.cache_action = (
            'matching maximum already complete'
            if float(metadata['tmax']) == target_tmax
            else 'resumed checkpoint to larger tmax'
        )
        return simulation

    @classmethod
    def run_cached(cls, data_dir='data', preferred_output=None,
                   overwrite=False, **parameters):
        """Reuse compatible C++ data, run only missing time, and persist it."""

        requested_tmax = parameters.get('tmax', 20.0)
        simulation = cls.from_data_folder(data_dir, **parameters)
        if (
            simulation.cache_action == 'matching maximum already complete'
            and simulation.time >= simulation.parameters['tmax']
        ):
            return (
                simulation, simulation.saved_rows,
                simulation.source_csv_path, simulation.source_json_path,
            )
        rows = simulation.run()
        if (
            preferred_output is not None
            and simulation.parameters['tmax'] == requested_tmax
        ):
            csv_path = Path(preferred_output)
        else:
            csv_path, _ = simulation.result_paths(data_dir)
        csv_path, json_path = simulation.write_results(csv_path, overwrite)
        cls._discard_superseded(simulation, csv_path, json_path)
        return simulation, rows, csv_path, json_path

    @staticmethod
    def _discard_superseded(simulation, csv_path, json_path):
        """Delete the shorter run this one resumed from, once it is replaced.

        Its rows are a strict prefix of the file just written, so keeping it
        only leaves two results for one set of parameters. Only the pair the
        run actually restored from is removed, and only after the new pair is
        safely on disk: an interrupted extension therefore loses nothing, and
        two concurrent jobs cannot delete each other's output.
        """

        written = {csv_path.resolve(), json_path.resolve()}
        for stale in (simulation.source_csv_path, simulation.source_json_path):
            if stale is not None and Path(stale).resolve() not in written:
                Path(stale).unlink(missing_ok=True)

    def diagnostics(self):
        """Return current C++ engine clocks, intervals, and counters."""

        return _engine.info(self._simulation)

    @property
    def time(self):
        return self.diagnostics()['time']

    @property
    def record_interval(self):
        return self.diagnostics()['record_interval']


class NNNHeap(NNN):
    """C++ engine keeping every process on one clock, ordered by a heap."""

    def __init__(self, *args, **kwargs):
        if 'dynamics' in kwargs and kwargs['dynamics'] != 'heap':
            raise ValueError('NNNHeap requires heap dynamics')
        kwargs['dynamics'] = 'heap'
        super().__init__(*args, **kwargs)


class NNNUniformized(NNN):
    """C++ baseline retaining blind proposals and blocked-event rejection."""

    def __init__(self, *args, **kwargs):
        if 'dynamics' in kwargs and kwargs['dynamics'] != 'uniformized':
            raise ValueError('NNNUniformized requires uniformized dynamics')
        kwargs['dynamics'] = 'uniformized'
        super().__init__(*args, **kwargs)


class NNNDirect(NNN):
    """C++ engine evolving each walker alone between fixed cloning instants."""

    def __init__(self, *args, **kwargs):
        if 'dynamics' in kwargs and kwargs['dynamics'] != 'direct':
            raise ValueError('NNNDirect requires direct dynamics')
        kwargs['dynamics'] = 'direct'
        super().__init__(*args, **kwargs)
