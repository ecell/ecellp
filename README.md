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

### Install as python library #NOT recommended
```
python setup.py install --prefix=${PREFIX}
```

### Install for virtualenv
Install virtualenv
```
pip install viratualenv virtualenvwrapper
```
Added the following in `.zshrc` or `.bashrc`
```bash
export WORKON_HOME=$HOME/.virtualenvs
if [ -f path_to_your_env/bin/virtualenvwrapper.sh ]; then
   source path_to_your_env/bin/virtualenvwrapper.sh
fi
```
Into virtualenv and add `${PYTHONPATH}` to this project directory
```
mkvirtualenv ecellp
workon ecellp
add2virtualenv "/home/soh.i/E-cell_Sprint/ecellp/src"
```

## Run test
```bash
python setup.py test
```

## Running
Sample code is `query_test.py` in `samples` directory.

```
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages #not necessary for virtualenv
python samples/query_test.py
```

### Configuration file
Add **relative path** from APP_ROOT to your data files.
```ini
[input_data]
sequence=/data/test.fa
```

### Query genome sequence and annotations
```python
from ecellp import session

db_config = session.DBConfig()

# pass absolute path to your conf.ini
# db_config = session.DBConfig(file=os.path.abspath('./conf.ini'))

# pass a dict containing paths to your data
# db_config = session.DBConfig(paths={'sequence':'/data/test.fa'})

query = session.QueryBuilder(db_config)

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
