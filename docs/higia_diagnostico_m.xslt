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
        <h1>Repositorio de Resultados</h1>

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

        <!-- Tabla 2: Características morfológicas -->
        <h2>Características Morfológicas de Plasmodium</h2>
        <table>
          <tr>
            <th>Tipo</th>
            <th>Glóbulo Rojo</th>
            <th>Trofozoito Joven</th>
            <th>Trofozoito Mediano</th>
            <th>Trofozoito Adulto</th>
            <th>Esquizonte</th>
            <th>Gametocito</th>
          </tr>
          <xsl:for-each select="determinaciones/caracteristicasMorfologicas/plasmodium">
            <tr>
              <td><xsl:value-of select="@tipo" /></td>
              <td><xsl:value-of select="globuloRojo" /></td>
              <td><xsl:value-of select="trofozoitoJoven" /></td>
              <td><xsl:value-of select="trofozoitoMediano" /></td>
              <td><xsl:value-of select="trofozoitoAdulto" /></td>
              <td><xsl:value-of select="esquizonte" /></td>
              <td><xsl:value-of select="gametocito" /></td>
            </tr>
          </xsl:for-each>
        </table>

        <!-- Tabla 3: Rangos de densidad parasitaria -->
        <h2>Rangos de Densidad Parasitaria</h2>
        <table>
          <tr>
            <th>Tipo</th>
            <th>Rango</th>
            <th>Por 100 leucocitos</th>
            <th>Por microlitro</th>
          </tr>
          <xsl:for-each select="determinaciones/rangosDensidadParasitaria/plasmodium">
            <xsl:variable name="tipo" select="@tipo"/>
            <xsl:for-each select="rango">
              <tr>
                <td><xsl:value-of select="$tipo" /></td>
                <td><xsl:value-of select="@densidad" /></td>
                <td><xsl:value-of select="por100leucocitos" /></td>
                <td><xsl:value-of select="pormicrolitro" /></td>
              </tr>
            </xsl:for-each>
          </xsl:for-each>
        </table>

      </body>
    </html>
  </xsl:template>

</xsl:stylesheet>




