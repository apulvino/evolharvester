##### setup.py
from setuptools import setup, find_packages

setup(
	name="evolharvester",
	version="0.1.0",
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	install_requires=["pandas"],
	entry_points={
		"console_scripts": [
			"evolharvester = evolharvester.cli:main",
		],
	},
)
