# Optica-HOE

Prototipo de un elemento óptico holográfico reflectivo para generar una imagen virtual a 2 m del observador.

## Ejecución

MATLAB R2018b o posterior. Colocar `parabrisas.stl` (solo el vidrio, unidades mm) y opcionalmente `logo.png` (figura clara sobre fondo oscuro) junto a `prototipo_HOE.m`. Abrir la carpeta en MATLAB y ejecutar `prototipo_HOE`. Sin esos archivos se usan geometría ilustrativa y texto VITRO de prueba. También se pueden indicar rutas: `prototipo_HOE('parabrisas.stl','logo.png')`.

## Resultados

La carpeta `resultados_HOE` contiene reconstruccion.png, integracion_CAD.png y datos_HOE.mat. Parámetros iniciales: 532 nm, HOE de 20 mm de lado, ojo a 600 mm, plano virtual a 1400 mm detrás del HOE, incidencia a 30 grados. Gerchberg-Saxton y propagación escalar; la portadora oblicua se trata analíticamente.

## Limitaciones

Versión inicial sin prueba de ejecución en MATLAB. El CAD sirve para posicionar un plano tangente; no hay trazado dentro del laminado ni compensación por curvatura. La distancia calculada con rayos ideales es una comprobación de referencia, no una medición independiente de profundidad del logo. Verificar escala, cara seleccionada y lado de la normal del STL. No representa una máscara fabricable ni eficiencia de un HOE volumétrico. No publicar el CAD del socio formador sin permiso.
