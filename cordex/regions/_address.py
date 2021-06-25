import requests
import tempfile


def VG2500(domain="lan"):
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
    url = "https://daten.gdz.bkg.bund.de/produkte/vg/vg2500/aktuell/vg2500_01-01.gk3.shape.zip"
    shp_file = "!vg2500_01-01.gk3.shape/vg2500/vg2500_{}.shp".format(domain)
    tmp = _download(url, ".zip")
    return "zip://" + tmp + shp_file


def _download(url, suffix=None):
    response = requests.get(url)
    f = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    f.write(response.content)
    f.close()
    return f.name
