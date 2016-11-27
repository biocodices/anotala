from setuptools import setup, find_packages


with open('requirements.txt') as f:
    requirements = f.read().split('\n')

with open('README.md') as f:
    long_description = f.read()

setup(
    name='anotamela',
    version='1.0',
    author='Juan Manuel Berros',
    author_email='juanma.berros@gmail.com',
    url='https://github.com/biocodices/anotamela',
    license='MIT',
    install_requires=requirements,
    packages=find_packages(),
)
