# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file is the basic package setup. It's currently fleshed out with a very
# minimal package list. If this goes into the PyPI repository, we need to add
# appropriate metadata, etc.

# 2021.02.25.: switched to setuptools @Elija, Feigl
from importlib.metadata import entry_points
from setuptools import setup, find_packages

setup(
    name="nanodesign",
    version="1.17",
    packages=find_packages(),
    include_package_data=True,
    install_requires=(
        "numpy",
        "biopython",
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Apache, Version 2.0",
        "Operating System :: OS Independent",
    ],
    entry_points="""
        [console_scripts]
        nanodesign=scripts.main:cli
    """,
)
