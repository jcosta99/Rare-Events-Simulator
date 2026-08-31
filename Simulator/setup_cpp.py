"""Build the C++ engine in the source tree.

Run this from the ``Simulator/`` folder that contains it: setuptools resolves
the package directory it copies the extension into against the working
directory, not against this file.

The standard build targets x86-64-v2. This lets the compiler translate the hot
``std::popcount`` calls into hardware POPCNT instructions on modern x86-64
laptops. See ``docs/CPU_REQUIREMENTS.md`` for the portable fallback.

Both intermediate products are written straight into ``build/``: the object
file and the staged extension, rather than the platform-tagged directory tree
setuptools uses by default. The importable copy still lands beside the Python
package, because the build is run with ``--inplace``.
"""

from pathlib import Path

from setuptools import Extension
from setuptools import setup
from setuptools.command.build_ext import build_ext


BUILD_DIR = Path(__file__).parent / 'build'


class FlatBuildExt(build_ext):
    """Keep the two build products directly inside ``build/``."""

    def initialize_options(self):
        super().initialize_options()
        self.build_temp = str(BUILD_DIR)
        self.build_lib = str(BUILD_DIR)

    def build_extensions(self):
        # Object files normally mirror the source tree, giving build/cpp/*.o.
        # Stripping the directory puts them straight in build/ -- and because
        # the same call decides where an existing object is looked for, an
        # unchanged source is still not recompiled.
        original = self.compiler.object_filenames

        def flat_object_filenames(sources, strip_dir=0, output_dir=''):
            return original(sources, strip_dir=1, output_dir=output_dir)

        self.compiler.object_filenames = flat_object_filenames
        super().build_extensions()

    def get_ext_fullpath(self, ext_name):
        # The staged extension normally mirrors the package tree, giving
        # build/Python_Cpp_Interface/_engine*.so.  Flatten it -- but only
        # inside build/, never the in-place copy, which belongs beside the
        # package.  This is the same path setuptools tests for freshness, so
        # flattening here rather than moving the file afterwards keeps
        # rebuilds incremental.
        path = Path(super().get_ext_fullpath(ext_name))
        if BUILD_DIR in path.parents:
            return str(BUILD_DIR / path.name)
        return str(path)

    def copy_extensions_to_source(self):
        # setuptools looks for the staged file under a package directory, so
        # the flattened layout needs its own copy step.  This is what puts the
        # importable extension beside Python_Cpp_Interface/__init__.py.
        build_py = self.get_finalized_command('build_py')
        for extension in self.extensions:
            fullname = self.get_ext_fullname(extension.name)
            filename = Path(self.get_ext_filename(fullname)).name
            package = '.'.join(fullname.split('.')[:-1])
            destination = Path(build_py.get_package_dir(package)) / filename
            self.copy_file(str(BUILD_DIR / filename), str(destination),
                           level=self.verbose)


setup(
    name='monitored-exclusion-cpp-engine',
    version='0.1.0',
    cmdclass={'build_ext': FlatBuildExt},
    ext_modules=[
        Extension(
            'Python_Cpp_Interface._engine',
            sources=['cpp/engine.cpp'],
            language='c++',
            extra_compile_args=['-O3', '-std=c++20', '-march=x86-64-v2'],
        ),
    ],
)
