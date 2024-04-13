#! /usr/bin/env python3

# Copyright 2018, University of New Hampshire

import sqlite3


class DataStore:
    """ Manages database connections and queries """
    _entries = dict()

    @staticmethod
    def destroy():
        """ Close all DB connections """
        for entry in DataStore._entries:
            entry.close()

    @staticmethod
    def get_entry(name):
        """ Get datastore entry """
        return DataStore._entries.get(name, None)

    def __init__(self, name, filename):
        self.queries = dict()
        self.indices = dict()
        self.connection = sqlite3.connect(filename)
        self.connection.isolation_level = None
        self.transaction = False

        # Create age table (if applicable)
        cursor = self.connection.cursor()
        cursor.execute("CREATE TABLE IF NOT EXISTS age (name text PRIMARY KEY, modified datetime)")
        self.connection.commit()

        # Add connection to DataStore
        DataStore._entries[name] = self

    def create_table(self, name, columns):
        """ Create a table given a list of (name, type, modifier) tuples """
        cursor = self.connection.cursor()
        create_ddl = ", ".join([" ".join(x) for x in columns])
        cursor.execute("CREATE TABLE IF NOT EXISTS {0} ({1})".format(name, create_ddl))

        if not self.transaction:
            self.connection.commit()

    def define_index(self, name, table, columns, unique):
        """ Define a column index """
        self.indices[name] = (table, columns, unique)

    def create_index(self, name):
        """ Create a column index """
        if name not in self.indices:
            return

        table, columns, unique = self.indices[name]

        cursor = self.connection.cursor()
        create_ddl = ", ".join(columns)
        unique_ddl = "UNIQUE" if unique else ""
        cursor.execute("CREATE {0} INDEX IF NOT EXISTS {1} ON {2}({3})".format(unique_ddl, name, table, create_ddl))

    def drop_index(self, name):
        """ Delete a column index """
        if name not in self.indices:
            return

        cursor = self.connection.cursor()
        cursor.execute("DROP INDEX IF EXISTS {0}".format(name))

    def insert_rows(self, name, data):
        """ Insert rows given a list of tuples """
        cursor = self.connection.cursor()
        param_tokens = ("?," * len(data[0]))[:-1]
        cursor.executemany("INSERT OR IGNORE INTO {0} VALUES ({1})".format(name, param_tokens), data)

        if not self.transaction:
            self.connection.commit()

    def delete_rows(self, name, query=None):
        """ Delete rows given a filter query """
        cursor = self.connection.cursor()
        where = " WHERE {0}".format(query) if query else ""
        cursor.execute("DELETE FROM {0}{1}".format(name, where))

        if not self.transaction:
            self.connection.commit()

    def define_query(self, name, query):
        """ Define a query """
        self.queries[name] = query

    def exec_query(self, name, parameters=""):
        """ Execute a query and retrieve the results if applicable """
        # Check for valid query
        if name not in self.queries:
            return None

        cursor = self.connection.cursor()
        result = cursor.execute(self.queries[name], parameters)

        if not self.transaction:
            self.connection.commit()

        return result

    def process_trans(self):
        """ Start/end transaction """
        cursor = self.connection.cursor()

        if not self.transaction:
            cursor.execute("BEGIN TRANSACTION")
        else:
            cursor.execute("END TRANSACTION")

        self.transaction = not self.transaction

    def update_age(self, name):
        """ Mark table age with current datetime """
        cursor = self.connection.cursor()
        cursor.execute("INSERT OR IGNORE INTO age VALUES(?, 0)", [name])
        cursor.execute("UPDATE age SET modified = CURRENT_TIMESTAMP WHERE name = ?", [name])
        self.connection.commit()

    def get_expired(self, name, days):
        """ Check if table is expired """
        cursor = self.connection.cursor()
        cursor.execute("SELECT JulianDay(CURRENT_TIMESTAMP) - JulianDay(modified) FROM age WHERE name = ?", [name])
        result = cursor.fetchone()

        if not result or result[0] > days:
            return True

        return False

    
