from setuptools import setup, find_packages


with open('requirements.txt') as f:
    requirements = f.read().split('\n')

with open('README.md') as f:
    long_description = f.read()

setup(
    name='anotala',
    version='1.0.1',
    author='Juan Manuel Berros',
    author_email='juanma.berros@gmail.com',
    url='https://github.com/biocodices/anotala',
    license='MIT',
    install_requires=requirements,
    packages=find_packages(),
)

