#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from setuptools import setup

# setup.
setup(name='brain_tumor_em',
      version='0.1',            # version is kind of arbitrary.
      description='Expectation-Maximisation algoirthm to segment '
        'brain tumors in MR scans',
      author='Esther Alberts',
      author_email='esther.alberts@tum.de',
      license='MIT',
      download_url='',
      packages=['em', 'utils','example'],
      install_requires=['numpy',
                        'SimpleITK'],

      # CLI
      scripts=['example/main.py'],
      )