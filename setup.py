from setuptools import setup, find_packages


setup(
    version='0.0.1',
    description='E-Cell P: Database project for whole cell simulation in E. coli',
    author='Soh Ishiguro',
    author_email='t10078si@sfc.keio.ac.jp',
    url='https://github.com/soh-i/DB4E-Cell-P/',
    name='ecellp',
    license='not yet',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=['biopython==1.62', 'sqlalchemy==0.8.2', 'pysqlite==2.6.3'],
    test_suite="ecellp.tests"
    )
