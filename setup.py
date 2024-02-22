import os
import setuptools

about = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'pagn', '__version__.py'), 'r') as _:
    exec(_.read(), about)

with open ('requirements.txt', 'r') as _:
    requires = [line.split()[0] for line in _]

setuptools.setup(
    long_description='See: '+about['__url__'],
    long_description_content_type='text/markdown',
    version=about['__version__'],
    url=about['__url__'],
    packages=setuptools.find_packages(),
    zip_safe=False
)
