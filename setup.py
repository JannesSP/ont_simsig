from setuptools import setup, find_packages

requirements = [
      'python',
      'h5py',
      'biopython',
      'matplotlib',
      'numpy',
      'scipy',
      'pandas',
      'seaborn',
      ]

setup(name='ont_simsig',
      version='1.0.0',
      description='Simulation of Oxford Nanopore Technologies RNA sequencing signals',
      author='Jannes Spangenberg',
      author_email='jannes.spangenberg@uni-jena.de',
      url='https://github.com/JannesSP/ont_simsig',
      license='GNU General Public License v3.0',
      packages=find_packages(),
      py_modules=['src'],
      package_data={
        'ont_simsig':['README.md','LICENSE'],
      },
      entry_points={
            'console_scripts':[
                  'ont_simsig=src.ont_simsig:main',
            ]
      },
      install_requires=requirements,
      keywords=['ont_simsig', 'ONT simulation', 'RNA simulation', 'ONT', 'signal simulation', 'raw data', 'Oxford Nanopore Technologies', 'MinION', 'Direct RNA Sequencing'],
     )
