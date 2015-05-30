#!/usr/bin/env python

from ConfigParser import SafeConfigParser
import logging as log
import os
import re
import sys
import warnings

sys.path.append("%s/lib/python" % os.environ['LAB_HOME'])
import util.cl
from util.collections import VarArgGenerator
from util.db import DBWrapper

def exec_sql(sql, config, sql_vars={}, ignore=[], dryrun=False, result_handler=None):
    stmts = map(lambda s: s.strip(), re.split(';', sql))
    # remove comments
    stmts = map(lambda s: re.sub('^#.*?\n+', '', s, flags=re.M), stmts)
    # remove empty statements
    stmts = filter(lambda s: len(s) > 0 and re.match('^[^\s\n]', s), stmts)
    # interpolate
    stmts = [s % sql_vars for s in stmts]
    
    if dryrun:
        log_stmts(stmts)
    else:
        exec_stmts(stmts, config, ignore, result_handler)

def log_stmts(stmts):
    for s in stmts: 
        log.info(s)

def exec_stmts(stmts, config, ignore, result_handler):
    with DBWrapper(config.pop('engine', 'mysql'), **config) as db:
        for s in stmts: 
            exec_stmt(s, db, ignore, result_handler)
            
def exec_stmt(stmt, db, ignore=[], result_handler=None):
    stmt_type = stmt[0:stmt.find(' ')]
    try:
        if stmt_type == 'SELECT' and not re.search('OUTFILE', stmt, re.I):
            results = db.selectall(stmt, ignore)
            if result_handler:
                result_handler(results)
        else:
            db.execute(stmt, ignore=ignore)
    except:
        log.error("Error executing statement %s" % stmt, exc_info=True)
        sys.exit()

if __name__ == '__main__':
    def add_arguments(parser):
        parser.add_argument('--db_config_file', type='readable_file', metavar="DIR", 
            default='%s/db/db.properties' % os.environ['LAB_HOME'], 
            help="File contining DB configuration parameters.")
        parser.add_argument('-d', '--database', metavar="NAME", default='mouse',
            help="Name of database to connect to.")
        parser.add_argument('-i', '--infile', action='dict', type='mapping', 
            default={}, metavar="NAME=PATH", 
            help="Define an input file or glob variable.")
        parser.add_argument('-o', '--outfile', action='dict', type='mapping', 
            default={}, metavar="NAME=PATH", 
            help="Define an output file variable.")
        parser.add_argument('-v', '--variable', action='dict', type='mapping', 
            metavar="NAME=VALUE", default={},
            help="Define a string variable to be replaced in SQL statements.")
        parser.add_argument('--dryrun', action='store_true', default=False,
            help="Print SQL to be executed, but don't execute it.")
        parser.add_argument('--ignore', type='str_list', action='extend', default=[], 
            help="Ignore: warnings, errors, skipped.")
        parser.add_argument('--overwrite', action='store_true', default=False,
            help="Overwrite exiting output files.")
        parser.add_argument('-s', '--sql', metavar="SQL", help="SQL to execute.")
        parser.add_argument('sql_file', type='readable_file', metavar="FILE", 
            nargs='?', help="File containing SQL to execute.")
    
    ns = util.cl.parse(add_arguments)
    
    config = SafeConfigParser()
    config.read(ns.db_config_file)
    if config.has_section(ns.database):
        config = dict(config.items(ns.database))
    else:
        log.error("Database '%s' not found in config file" % ns.database)
        sys.exit()
    
    sql_vars = ns.variable or {}
    
    for name,path in ns.infile:
        try:
            sql_vars[name] = util.io.check_path(util.io.resolve_path(path), 'f', 'r')
        except IOError:
            log.error(exc_info=True)
            sys.exit()

    for name,path in ns.outfile:
        path = os.path.expanduser(path)
        try:
            path = util.io.abspath(path)
            if os.path.exists(path):
                if ns.overwrite:
                    os.remove(path)
                else:
                    log.error("File already exists %s" % path)
            else:
                parent = os.path.dirname(path)
                if not os.path.exists(parent):
                    os.makedirs(parent)
            sql_vars[name] = path
        except IOError:
            log.error(exc_info=True)
            sys.exit()
        
    sql = ns.sql or util.io.safe_read_file(ns.sql_file)
    exec_sql(sql, config, sql_vars, set(ns.ignore), ns.dryrun)