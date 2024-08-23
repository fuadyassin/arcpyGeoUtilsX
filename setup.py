from setuptools import setup, find_packages

setup(
    name='arcpyGeoUtilsX',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # Add other dependencies here
    ],
    include_package_data=True,
    description='A Python package for GIS processing using arcpy.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/fuadyassin/arcpyGeoUtilsX',
    author='Fuad Yassin',
    author_email='fuad.yassin@usask.ca',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
