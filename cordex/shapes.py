from . import data

"""convert to WGS84 latitude-longitude projection"""
WGS84 = "EPSG:4326"


def get_geodataframe(shapefile, to_crs=WGS84, **kwargs):
    import geopandas as gpd

    shp = gpd.read_file(shapefile)
    if to_crs is not None:
        shp = shp.to_crs(to_crs)
    return shp


def get_regionmask(geodataframe, **kwargs):
    import regionmask

    return regionmask.from_geopandas(geodataframe, **kwargs)


class VG2500:
    """VG2500 Deutschland Verwaltungsgrenzen

    ADE Administrative Ebene
        Werteübersicht: 1 = Staat 2 = Land 3 = Regierungsbezirk 4 = Kreis
    ARS Amtlicher Regionalschlüssel (bisher Attribut RS)
        Bei diesem Schlüssel handelt es sich um den statistischen Schlüssel. Der Schlüssel ist hierarchisch
        strukturiert und spiegelt die in der Bundesrepublik Deutschland bestehenden Verwaltungsebenen wider.
        Der ARS gliedert sich wie folgt:
        1. – 2. Stelle = Kennzahl des Landes
        3. Stelle = Kennzahl des Regierungsbezirks
        4. – 5. Stelle = Kennzahl des Kreises
        6. – 9. Stelle = Kennzahl der Verwaltungsgemeinschaft
        10. – 12. Stelle = Kennzahl der Gemeinde
    GEN Geografischer Name
    ARS_0 aufgefüllter Amtlicher Regionalschlüssel (bisher Attribut RS_0)
        grundsätzlich 12-stelliger ARS (mit Nullen rechtsseitig aufgefüllt)
    """

    @staticmethod
    def _filename(domain):
        # url = "https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/vg2500_01-01.gk3.shape.zip"
        shp_file = "!vg2500_01-01.gk3.shape/vg2500/vg2500_{}.shp".format(domain)
        fname = data.fetch_vg2500()
        return "zip://" + fname + shp_file

    @classmethod
    def geodata(cls, domain="lan"):
        """Returns a GeoDataFrame object."""
        url = cls._filename(domain)
        geodata = get_geodataframe(url)
        geodata["name"] = geodata["ARS"] + "_" + geodata["GEN"]
        return geodata

    @classmethod
    def regionmask(cls, domain="lan"):
        """Returns a maks."""
        return get_regionmask(cls.geodata(domain), names="name", abbrevs="_from_name")
