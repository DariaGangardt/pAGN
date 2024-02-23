import os
import setuptools

about = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'pagn', '__version__.py'), 'r') as _:
    exec(_.read(), about)

setuptools.setup(
    version=about['__version__'],
)
