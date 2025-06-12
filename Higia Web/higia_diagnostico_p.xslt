<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="html" encoding="UTF-8" indent="yes" />

  <xsl:template match="/datos">

    <html>
      <head>
        <title>Resultados Diagnóstico Higia</title>
		<link rel="stylesheet" href="higiausuario.css" type="text/css"/>
      </head>

      <body>

        <!-- Tabla 1: Diagnósticos -->
        <h2>Diagnósticos de Pacientes</h2>
        <table style="width: 95%; font-size: 13px; table-layout: auto; white-space: nowrap;">
          <tr>
            <th>ID</th>
            <th>Nombre</th>
            <th>Apellido</th>
            <th>Sexo</th>
            <th>Nacimiento</th>
            <th>Fecha</th>
            <th>Imagen</th>
            <th>Resultado</th>
            <th>Especie</th>
            <th>Densidad</th>
          </tr>
          <xsl:for-each select="diagnosticos/diagnostico">
            <tr>
              <td><xsl:value-of select="@id" /></td>
              <td><xsl:value-of select="paciente/nombre" /></td>
              <td><xsl:value-of select="paciente/apellido" /></td>
              <td><xsl:value-of select="paciente/sexo" /></td>
              <td><xsl:value-of select="paciente/nacimiento" /></td>
              <td><xsl:value-of select="fecha" /></td>
				<td>
				  <img 
					alt="Frotis" 
					width="100" 
					style="cursor: zoom-in;">
					<xsl:attribute name="src">
					  <xsl:value-of select="datosimagen/imagen"/>
					</xsl:attribute>
					<xsl:attribute name="onclick">
					  <xsl:text>parent.ampliarImagen(this.src)</xsl:text>
					</xsl:attribute>
				  </img>
				</td>
              <td><xsl:value-of select="resultado/resultadotest" /></td>
              <td><xsl:value-of select="resultado/especiePlasmodium" /></td>
              <td><xsl:value-of select="resultado/densidadParasitaria" /></td>
            </tr>
          </xsl:for-each>
        </table>

      </body>
    </html>
  </xsl:template>

</xsl:stylesheet>
