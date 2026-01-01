from setuptools import setup, find_packages

setup(
    name='multiscale2',
    version='0.1.0',
    packages=find_packages(),
    description='A multi-stage workflow for simulating protein condensates.',
    author='Xiaojing Tian',
    author_email='tianxj15@tsinghua.org',
    entry_points={
        'console_scripts': [
            'ms2=multiscale2.cli:main',
        ],
    },
    install_requires=[
        'click>=8.0',
        'pyyaml>=6.0',
    ],
)

