# Optica-HOE

Prototipo numérico de un elemento óptico holográfico reflectivo para reconstruir el logo de Vitro a 2 m del observador, con geometría CAD del parabrisas.

## Versión final — ejecución

Usar MATLAB Online o una instalación de MATLAB compatible con `stlread` (devuelve una triangulación), `xline` y `savefig` (se recomienda R2018b o posterior).

Colocar estos tres archivos en la misma carpeta:

- `HOE_Final.m` — programa autónomo.
- `parabrisas.stl` — CAD del vidrio, coordenadas en mm.
- `logo.png` — el mismo logo claro sobre fondo oscuro de las pruebas.

Ejecutar:

```matlab
HOE_Final
```

No requiere ejecutar las revisiones anteriores ni conservar sus archivos MAT. El CAD se ajusta de nuevo en cada ejecución. También acepta rutas mediante `HOE_Final(cadFile,logoFile)`.

## Parámetros definitivos

532 nm, incidencia 30°, HOE 20 × 20 mm, pupila de 5 mm, ojo a 600 mm del HOE, plano virtual 1400 mm detrás del HOE (2 m desde el ojo), logo de 50 mm de ancho, malla 2048 × 2048, semilla 7, 40 iteraciones convencionales y 40 de optimización con pupila.

## Estructura del código

1. Lectura y ajuste local del CAD.
2. Parámetros y preparación del logo.
3. Propagación de Fresnel.
4. Diseño iterativo de fase: Gerchberg-Saxton y optimización con pupila.
5. Barrido axial.
6. Compensación ideal de fase por curvatura.
7. Guardado de resultados.
8. Figuras ligeras sin transparencias ni texturas 3D.

## Archivos generados

En `resultados_Final/`:

- `datos_HOE.mat` (se guarda antes de generar gráficos).
- `reconstruccion.png`, `profundidad.png`, `curvatura.png`.
- `prototipo_3D.png`, `prototipo_3D.fig`.
- `geometria_3D.png`, `geometria_3D.fig`, `logo_observado.png`.

## Comprobación

La R3 de 50 mm ejecutada en MATLAB Online mostró captura de 37.68 %, similitud de amplitud de 0.901 y máximo axial a 2.00 m. Son valores de referencia redondeados, no resultados ya verificados de RFinal. Comparar la salida nueva con ellos usando los mismos archivos. RFinal conserva el núcleo óptico de R3 y el ajuste geométrico de R2, pero todavía requiere una ejecución de comprobación en MATLAB.

## Alcance y limitaciones

Modelo escalar idealizado de fase reflectiva en aire. No calcula refracción multicapa, eficiencia de difracción volumétrica ni fabricabilidad. La compensación de curvatura es una cancelación algebraica ideal, no validación independiente. Si el paso de fase supera pi, se omite la predicción cuantitativa del caso sin compensar. La malla STL se simplifica solamente para dibujar; comprobar orientación, unidades y que los puntos del ajuste pertenecen al vidrio. La captura es potencia relativa del campo ideal que entra en la pupila, no eficiencia de un dispositivo real. La gráfica 3D es geométrica; la imagen del logo procede de la reconstrucción numérica, con gamma de presentación. No publicar el CAD del socio formador sin permiso.

## Historial de desarrollo

Las siguientes instrucciones son históricas y no son necesarias para ejecutar la versión final. Se conservan todos los scripts anteriores.

Prototipo de un elemento óptico holográfico reflectivo para generar una imagen virtual a 2 m del observador.

## Ejecución

MATLAB R2018b o posterior. Colocar `parabrisas.stl` (solo el vidrio, unidades mm) y opcionalmente `logo.png` (figura clara sobre fondo oscuro) junto a `prototipo_HOE.m`. Abrir la carpeta en MATLAB y ejecutar `prototipo_HOE`. Sin esos archivos se usan geometría ilustrativa y texto VITRO de prueba. También se pueden indicar rutas: `prototipo_HOE('parabrisas.stl','logo.png')`.

## Resultados

La carpeta `resultados_HOE` contiene reconstruccion.png, integracion_CAD.png y datos_HOE.mat. Parámetros iniciales: 532 nm, HOE de 20 mm de lado, ojo a 600 mm, plano virtual a 1400 mm detrás del HOE, incidencia a 30 grados. Gerchberg-Saxton y propagación escalar; la portadora oblicua se trata analíticamente.

## Limitaciones

Versión inicial sin prueba de ejecución en MATLAB. El CAD sirve para posicionar un plano tangente; no hay trazado dentro del laminado ni compensación por curvatura. La distancia calculada con rayos ideales es una comprobación de referencia, no una medición independiente de profundidad del logo. Verificar escala, cara seleccionada y lado de la normal del STL. No representa una máscara fabricable ni eficiencia de un HOE volumétrico. No publicar el CAD del socio formador sin permiso.

## R1 — Optimización con pupila

Ejecutar `prototipo_HOE_R1` con los mismos archivos locales. Conserva R0. Compara Gerchberg-Saxton reproducido con el mismo operador de propagación frente a optimización de fase con pupila y regularización de apertura completa (peso 0.15). Guarda `resultados_HOE_R1/` con reconstruccion.png, optimizacion.png, integracion_CAD.png y datos_HOE.mat. El barrido axial de 1.4 a 2.6 m compara similitud de amplitud en coordenadas físicas; su máximo es una estimación numérica, no una medición experimental. El CAD sigue siendo tangente, sin corrección de curvatura ni refracción. R1 pendiente de ejecución en MATLAB.

## R2 — Curvatura local (modelo simplificado)

Ejecutar `prototipo_HOE_R2` con el mismo STL y PNG y el archivo `resultados_HOE_R1/datos_HOE.mat`. Si no existe, ejecuta R1 automáticamente. Si cambian el STL o el logo, regenerar R1 antes. R2 ajusta una superficie cuadrática a puntos locales del STL y aplica una pantalla de fase reflectiva en aire: delta = -k(1+cos(theta))*h. Compara plano, curvo sin compensar y compensación ideal. La cancelación es algebraica por construcción, no una validación física independiente. No modela refracción, capas, variación de eficiencia ni fabricación. El STL debe tener suficientes puntos locales de una misma cara; si falla el ajuste, refinar la exportación. Si el paso de fase supera pi, el caso sin compensar está submuestreado y no es cuantitativo. Revisar RMS, selección de cara y escala. Resultados en `resultados_HOE_R2/`: comparacion.png, curvatura.png, prototipo.png y datos_HOE.mat. R2 pendiente de ejecución en MATLAB con el STL real.

## R2.1 — Presentación del prototipo

Ejecutar `prototipo_HOE_R21` después de R2, con el mismo STL en milímetros. No recalcula ni modifica la física. Usa la intensidad reconstruida desde la pupila guardada en R2 como textura transparente en el plano virtual. Vista global ampliada 25 veces (solo visualización); vista desde el ojo con escala física y zoom angular. Se puede cambiar la ampliación con `prototipo_HOE_R21('parabrisas.stl',1)`. Guarda figuras FIG y PNG en `resultados_HOE_R21/`. El render no es una simulación retinal ni una validación independiente. Advierte si la fase curva no compensada está submuestreada. Pendiente de prueba en MATLAB Online.

## R3 — Prototipo de tamaño configurable

Ejecutar `prototipo_HOE_R3` (50 mm nominal) o `prototipo_HOE_R3(80)` (prueba de campo extendido, más memoria). Requiere el mismo `parabrisas.stl`, `logo.png` y `resultados_HOE_R2/datos_HOE.mat` del mismo CAD para reutilizar el ajuste local. Conserva 532 nm, 30°, HOE 20 × 20 mm, pupila 5 mm y distancias 0.6/1.4 m. Recorta márgenes negros del PNG para definir el ancho físico del contenido. Propaga hasta la pupila con convolución lineal de Fresnel, y mantiene el diseño iterativo con adjunto. El ancho geométrico sin viñeteo es aproximadamente 55 mm; 80 mm es una prueba, no un resultado garantizado. Salidas separadas en `resultados_HOE_R3_50mm/` etc. La visualización 3D usa la intensidad calculada a escala física; el zoom y gamma son solo visuales. La curvatura compensada es ideal por construcción; no hay refracción multicapa ni modelo volumétrico. Se omite la comparación cuantitativa sin compensar si el paso de fase excede pi. Pendiente de ejecución y convergencia en MATLAB Online.
