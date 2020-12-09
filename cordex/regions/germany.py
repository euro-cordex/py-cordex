
from . import _address as addr
from . import mask

"""convert to WGS84 latitude-longitude projection"""
WGS84 = "EPSG:4326"



class VG2500():
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
    def geodata(domain='lan'):
        """
        """
        url = addr.VG2500(domain)
        geodata = mask.get_geodata(url)
        geodata['name'] = geodata["ARS"] + '_' + geodata['GEN']
        return geodata

    @classmethod
    def regionmask(cls, domain='lan'):
        return mask.get_regionmask(cls.geodata(domain), names="name", abbrevs='_from_name')

