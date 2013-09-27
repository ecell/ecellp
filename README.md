E-Cell4 P: Database project for whole cell simulation in _E. coli_
===========
## Description 
Providing generating database class and query/filter DB interfaces for running simulation model

## Installation

### Install dependencies
```
sudo apt-get install libsqlite3-dev

pip install --install-option="--prefix=${PREFIX}" -r requeirements.txt
```

### Configuration file
Add **absolute** project root path into `conf.ini`.
```ini
[root]
APP_ROOT=/Users/yukke/dev/DB4E-Cell-P
```

### Install the library
```
python setup.py install --prefix=${PREFIX}
```

## Running
Sample code is `query_test.py` in `samples` directory.

```
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages python samples/query_test.py
```

### Query genome sequence and annotations
```python
from db4ecellp import db4ecellp
query = db4ecellp.QueryBuilder('./conf.ini') # PATH TO conf.ini

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
