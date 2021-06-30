from ._resources import read_cordex_tables, region_tables


class domains_cls():

    @property
    def tables(self):
        return read_cordex_tables()


domains = domains_cls()
