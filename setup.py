from setuptools import setup

setup(
    install_requires=[
        "numpy >= 1.13",
        "pandas >= 0.21.0",
        "tqdm >= 4.19",
        "json >= 1.35",
        "matplotlib >= 2.11"
    ],
    name='bashthebug',
    version='0.1.0',
    url='https://github.com/philipwfowler/bashthebug',
    author='Philip W Fowler',
    packages=['bashthebug'],
    license='MIT',
    scripts=['bin/analyse-bashthebug-classifications.py'],
    long_description=open('README.md').read(),
)
