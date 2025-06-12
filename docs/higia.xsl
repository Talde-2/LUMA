<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
    
    <xsl:output method="html" encoding="UTF-8" indent="yes"/>
    
    <xsl:template match="/">
        <html>
        <head>
            <title>Resultados Diagnósticos</title>
            <style>
                body { font-family: Arial; background-color: #f4faff; padding: 20px; }
                h2 { color: #004aad; }
                table { border-collapse: collapse; width: 100%; margin-bottom: 30px; }
                th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }
                th { background-color: #dbe9ff; color: #004aad; }
            </style>
        </head>
        <body>
            <h2>Resultados de Diagnósticos</h2>
            <table>
                <tr>
                    <th>ID</th>
                    <th>Fecha</th>
                    <th>Resultado</th>
                    <th>Especie</th>
                    <th>Densidad Parasitaria</th>
                </tr>
                <xsl:for-each select="diagnosticos/diagnostico">
                    <tr>
                        <td><xsl:value-of select="@id"/></td>
                        <td><xsl:value-of select="fecha"/></td>
                        <td><xsl:value-of select="resultado/positivoNegativo"/></td>
                        <td><xsl:value-of select="resultado/especiePlasmodium"/></td>
                        <td><xsl:value-of select="resultado/densidadParasitaria"/></td>
                    </tr>
                </xsl:for-each>
            </table>

            <h2>Características Morfológicas</h2>
            <table>
                <tr>
                    <th>Especie</th>
                    <th>Glóbulo Rojo</th>
                    <th>Fase Anillo</th>
                    <th>Trofozoíto Tardío</th>
                    <th>Esquizonte Maduro</th>
                    <th>Gametocitos</th>
                </tr>
                <xsl:for-each select="diagnosticos/determinaciones/caracteristicasMorfologicas/especie">
                    <tr>
                        <td><xsl:value-of select="@nombre"/></td>
                        <td><xsl:value-of select="globuloRojo"/></td>
                        <td><xsl:value-of select="faseAnillo"/></td>
                        <td><xsl:value-of select="trofozoitoTardio"/></td>
                        <td><xsl:value-of select="esquizonteMaduro"/></td>
                        <td><xsl:value-of select="gametocitos"/></td>
                    </tr>
                </xsl:for-each>
            </table>

            <h2>Rangos de Densidad Parasitaria</h2>
            <table>
                <tr>
                    <th>Nivel</th>
                    <th>P. falciparum (por 100 leucocitos)</th>
                    <th>P. falciparum (por μL)</th>
                    <th>P. vivax (por 100 leucocitos)</th>
                    <th>P. vivax (por μL)</th>
                </tr>
                <xsl:for-each select="diagnosticos/determinaciones/rangosDensidadParasitaria/rango">
                    <tr>
                        <td><xsl:value-of select="@nivel"/></td>
                        <td><xsl:value-of select="pFalciparum/por100Leucocitos"/></td>
                        <td><xsl:value-of select="pFalciparum/porMicrolitro"/></td>
                        <td><xsl:value-of select="pVivax/por100Leucocitos"/></td>
                        <td><xsl:value-of select="pVivax/porMicrolitro"/></td>
                    </tr>
                </xsl:for-each>
            </table>
        </body>
        </html>
    </xsl:template>
</xsl:stylesheet>

