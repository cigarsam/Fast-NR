from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='FastNR',
    version='1.0',
    description='FastNR is a software which applied on STARR-seq data to identify the negative regulatory elements, like silencer or insulator.',
    long_description=readme,
    author='Na He',
    author_email='11849492@mail.sustech.educ.com',
    url='https://github.com/Na-He/FastNR',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    entry_points = {
          'console_scripts': [
              'FastNR = src.FastNR:main'
          ]
          }
)

