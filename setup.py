#!/usr/bin/env python

from setuptools import setup, find_packages

with open("requirements/install.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="mab-libs",
    version='0.0.0',
    description="mab-libs is a library for generating synthetic antibody libraries in-silico",
    author="William Steele",
    url="git@github.com:wjs20/repo-name.git",
    packages=find_packages(exclude=["tests", "tests.*"]),
    install_requires=requirements,
    include_package_data=True,
)