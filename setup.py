from setuptools import setup, find_packages
from os import path
import re

package_name = "capital"
root_dir = path.abspath(path.dirname(__file__))

with open(path.join(root_dir, "README.md"), "r") as fh:
    long_description = fh.read()


with open(path.join(root_dir, package_name, '__init__.py')) as f:
    init_text = f.read()
    version = re.search(
        r'__version__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    license = re.search(
        r'__license__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    author = re.search(
        r'__author__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    author_email = re.search(
        r'__author_email__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)
    url = re.search(r'__url__\s*=\s*[\'\"](.+?)[\'\"]', init_text).group(1)


setup(
    name=package_name,
    packages=find_packages(),
    version=version,
    license=license,
    install_requires=[
        name.rstrip() for name in
        open(path.join(root_dir, 'requirements.txt')).readlines()
        ],
    author=author,
    author_email=author_email,
    url=url,
    description='Single-Cell Analysis, comparing pseudotime trajectories \
                with tree alignment',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: BSD License',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
