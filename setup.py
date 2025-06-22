from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from setuptools import find_packages
from setuptools import setup

description="Analyzing anharmonicity of traps for charged particles"

setup(
    name="anharm_analysis",
    version="1.0.0",
    description="Library for analyzing charged-particle trap anharmonicity",
    long_description=description,
    author="Electron",
    keywords=["Quantum Information Processing", "Numerical Simulation"],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "future",
        "numpy",
        "scipy",
        "pandas",
        "setuptools",
        "matplotlib",
        'tqdm',
        # 'optuna',
        # 'plotly',
        'scikit-learn',
        #'sphericart',
        'pyvista',
        'xarray',
        'cvxpy',
        'sympy',
        'juliacall'
    ],
    package_data = {
        "anharm_analysis": ["*.py"]
    },
    setup_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3.8",
    ],
    scripts=[
    ],
)
