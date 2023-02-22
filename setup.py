from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
    'h5py',
    'biopython',
    'matplotlib',
    'numpy',
    'pandas',
    'seaborn',
]

setup(
    name='ont_simsig',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Simulation of Oxford Nanopore Technologies RNA sequencing signals",
    license="MIT",
    author="Jannes Spangenberg",
    author_email='jannes.spangenberg@uni-jena.de',
    url='https://github.com/JannesSP/ont_simsig',
    packages=['ont_simsig'],
    entry_points={
        'console_scripts': [
            'ont_simsig=ont_simsig.cli:main'
        ]
    },
    install_requires=requirements,
    keywords='ont_simsig',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
