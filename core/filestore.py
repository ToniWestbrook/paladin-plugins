#! /usr/bin/env python3

# Copyright 2018, University of New Hampshire

import gzip
import os
import shutil
import tarfile
import tempfile
import time
from urllib import request


class FileStore:
    """ Manage temporary and cache files """
    FTYPE_TEMP, FTYPE_CACHE, FTYPE_OUTPUT, FTYPE_USER = range(4)
    FOPT_NORMAL, FOPT_GZIP, FOPT_GZIP_DECOMPRESS, FOPT_TAR = range(4)

    # Entries and system generated temporary path
    _entries = dict()
    _output_base = ""
    _cache_path = ""
    _temp_path = ""
    _output_path = ""
    _expire_age = 0

#    def _populate():
#        """ Add all entries to the file store """
#        FileStore("input", "input", "input.fq", None, FileStore.FTYPE_TEMP, FileStore.FOPT_NORMAL)

    @staticmethod
    def init(output_base, cache_path, output_path, expire_age):
        """ Initialize the file store  """
        FileStore._output_base = os.path.expanduser(output_base)
        FileStore._cache_path = os.path.expanduser(cache_path)
        FileStore._temp_path = tempfile.mkdtemp(prefix=output_base)
        FileStore._output_path = os.path.expanduser(output_path)
        FileStore._expire_age = expire_age

        # Create normal directories
        if not os.path.exists(FileStore._output_path):
            os.makedirs(FileStore._output_path)
        if not os.path.exists(FileStore._cache_path):
            os.makedirs(FileStore._cache_path)

    @staticmethod
    def destroy():
        """ Destroy the file store """
        # Close all open files
        for group in FileStore._entries.values():
            for entry in group.values():
                if entry.handle:
                    entry.handle.close()

        # Delete temporary directorty
        shutil.rmtree(FileStore._temp_path)

    @staticmethod
    def get_entry(name, group=None):
        """ Get a file entry from the store """
        if group:
            return FileStore._entries[group][name]
        else:
            return FileStore._entries[name][name]

    @staticmethod
    def get_group(group):
        """ Get all entries for a group """
        return FileStore._entries[group].values()

    def __init__(self, group, fid, name, url, ftype, options):
        self.fid = fid
        self.url = url
        self.ftype = ftype
        self.options = options
        self.handle = None

        # Set full path
        self.set_path(name, ftype)

        # Add current entry to file store
        FileStore._entries.setdefault(group, dict())
        FileStore._entries[group][fid] = self

    def get_path(self):
        if self.options != FileStore.FOPT_GZIP_DECOMPRESS:
            return self.path

        # Check if file already decompressed (partial/interrupted possible)
        decompressed = self.path[:-3]
        if os.path.exists(decompressed):
            return decompressed
        else:
            return self.path

    def set_path(self, name, ftype):
        """ Calculate and set the full path associated with the file """
        if ftype == FileStore.FTYPE_TEMP:
            self.path = os.path.join(FileStore._temp_path, name)
        if ftype == FileStore.FTYPE_CACHE:
            self.path = os.path.join(FileStore._cache_path, name)
        if ftype == FileStore.FTYPE_OUTPUT:
            self.path = os.path.join(FileStore._output_path, name)
        if ftype == FileStore.FTYPE_USER:
            self.path = name

    def get_handle(self, mode=None):
        """ Return file handle, and (re)open if mode specified """
        if mode:
            # Close if previously open
            if self.handle:
                self.handle.close()

            if self.options == FileStore.FOPT_GZIP:
                self.handle = gzip.open(self.get_path(), mode)
            else:
                self.handle = open(self.get_path(), mode)

        return self.handle

    def prepare(self):
        """ Download and/or decompress the file """
        # Download file if URL present
        if self.url:
            request.urlretrieve(self.url, self.path)

        # Decompress if gzip marked for decompression
        if self.options == FileStore.FOPT_GZIP_DECOMPRESS:
            with gzip.open(self.path, "rb") as handle_in:
                with open(self.path[:-3], "wb") as handle_out:
                    handle_out.write(handle_in.read())

            # Remove old file, update path to decompressed file
            os.remove(self.path)
            self.path = self.path[:-3]

        # Extract if an archive
        if self.options == FileStore.FOPT_TAR:
            handle = tarfile.open(self.path, "rt")
            handle.extractall(os.path.dirname(self.path))

    def exists(self):
        """ Check if file exists """
        return os.path.exists(self.get_path())

    def expired(self):
        """ Check if file is expired """
        mtime = os.path.getmtime(self.get_path())
        age = time.time() - mtime
        return (age / 86400) > FileStore.EXPIRE_AGE
