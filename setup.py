from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='Fast-NR',
    version='1.0',
    description='Fast-NR is a software which applied on STARR-seq data to identify the negative regulatory elements, like silencer or insulator.',
    long_description=readme,
    author='Na He',
    author_email='11849492@mail.sustech.educ.com',
    url='https://github.com/Na-He/Fast-NR',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    entry_points = {
          'console_scripts': [
              'Fast-NR = src.FastNR:main'
          ]
          }
)

