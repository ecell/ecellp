E-Cell P: Database project for whole cell simulation in _E. coli_
===========
## Description 
Providing generating database class and query/filter DB interfaces for running simulation model

## Installation

### Install dependencies
```
sudo apt-get install libsqlite3-dev
pip install --install-option="--prefix=${PREFIX}" -r requeirements.txt # if non-root user
```

### Configuration file
Add **absolute** project root path into `conf.ini`.
```ini
[root]
APP_ROOT=/Users/yukke/dev/DB4E-Cell-P
```

### Install the library #NOT recommended
```
python setup.py test # Run test
python setup.py install --prefix=${PREFIX} 
```

### Install for virtualenv
Install virtualenv
```
pip install viratualenv virtualenvwrapper
```
Added the following in `.zshrc` or `.bashrc`
```shell
export WORKON_HOME=$HOME/.virtualenvs
if [ -f $HOME/.pythonbrew/pythons/Python-2.7.5/bin/virtualenvwrapper.sh ]; then
   source $HOME/.pythonbrew/pythons/Python-2.7.5/bin/virtualenvwrapper.sh
fi
```
Into virtualenv and add ${PYTHONPATH} to this project directory
```
mkvirtualenv ecellp
workon ecellp
add2virtualenv "/home/soh.i/E-cell_Sprint/ecellp/src"
```

## Running
Sample code is `query_test.py` in `samples` directory.

```
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages
python samples/query_test.py
```

### Query genome sequence and annotations
```python
from ecellp import session
query = session.QueryBuilder('./conf.ini') # PATH TO conf.ini

print query.count_stored_records()
#=> 4145

print query.collect_cds_records()
#=> return All CDS annotations ans its sequences

print query.find_by_name('thrL')
#=> <Species('thrL','1','189','255','CDS','ATGAAACGCATTACCACCACCAT...')>

print query.collect_annotations_filter_by_strand(-1)
#=> return annotations and sequence on complement strand

print query.include_gene_in_region(21, 2012)
#=> [u'thrL', u'thrA']
```
