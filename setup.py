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
    author=about['__author__'],
    author_email=about['__author_email__'],
    url=about['__url__'],
    project_urls={
        'Source Code': about['__url__'],
        'Documentation': 'https://DariaGangardt.github.io/pAGN'
    },
    packages=setuptools.find_packages(),
    package_dir = {"":"."},
    include_package_data=True,
    install_requires=requires,
    package_data= {'': ['./opacity_tables/*']},
    zip_safe=False
)
