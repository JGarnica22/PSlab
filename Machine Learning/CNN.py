
# Implementación de una red neuronal convolucional para clasificar radiografías en normales o con derrame.

# Cargamos librerias
import numpy as np
import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.datasets import mnist
from tensorflow.keras.utils import to_categorical, plot_model, image_dataset_from_directory
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv2D, Flatten, MaxPooling2D, Dropout, Rescaling
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
import PIL
import pathlib

# Vemos algunas imágenes de radiografías normales y con derrame. Primero debemos desempaquetar (unzip) los ficheros y luego mostrarlos. Ponemos las dos carpetas descomprimidas en una carpeta común para después montar el dataset.

# Convocamos al terminal para unzip los files
!unzip -q normal.zip
!unzip -q effusion.zip
!ls
!mkdir radiografia
!mv normal effusion radiografia/

# Al cargar las imágenes hacemos ya el reshape a 64x64x1 y creamos los datasets de training y testing.

# cargamos imágenes, las preprocesamos a un tamaño de 64x64x1 y creamos los sets de training y testing
image_size = (64, 64)
batch_size = 100

train_orig, test_orig = tf.keras.preprocessing.image_dataset_from_directory(
    "radiografia",
    validation_split=100/700, # dividimos conjunto en 600 para train y 100 para test
    subset="both",
    seed=123,
    image_size=image_size,
    batch_size=batch_size,
    color_mode="grayscale" # indicamos color gris para que sólo haya un canal
)


# visualización de 12 imágenes normales y con derrame del set de training:
plt.figure(figsize=(10, 10))
for images, labels in train_orig.take(1):
    for i in range(12):
        ax = plt.subplot(3, 4, i + 1)
        plt.imshow(images[i].numpy().astype("uint8"), cmap="gray")
        plt.title(int(labels[i]))
        plt.axis("off")


# comprobamos el formato creado
for images, labels in train_orig.take(1):
    print("Shape of the original train data: {}".format(images.shape))



# Normalizamos mediante min-max los datsets, de training y testing.


# valores originales:
for x, y in train_orig.take(1):
    min = int(tf.reduce_min(x))
    max = int(tf.reduce_max(x))
    print("Los valores mínimos y máximos antes de normalizar son {} y {}.".format(min, max))



#normalizamos mediante transformación min-max
normalization_layer = Rescaling((min+1)/max) # indicamos cómo será el rescaling

norm_train = train_orig.map(lambda x, y: (normalization_layer(x), y))
image_batch, labels_batch = next(iter(norm_train))

for x, y in norm_train.take(1):
    min = int(tf.reduce_min(x))
    max = int(tf.reduce_max(x))
    print("Los valores mínimos y máximos después de normalizar son {} y {}.".format(min, max))


# repetimos proceso con el set de testing
# valores originales:
for x, y in test_orig.take(1):
    min = int(tf.reduce_min(x))
    max = int(tf.reduce_max(x))
    print("Los valores mínimos y máximos antes de normalizar son {} y {}.".format(min, max))

normalization_layer = Rescaling((min+1)/max) # indicamos cómo será el rescaling
norm_test = test_orig.map(lambda x, y: (normalization_layer(x), y))
image_batch, labels_batch = next(iter(norm_test))

for x, y in norm_test.take(1):
    min = int(tf.reduce_min(x))
    max = int(tf.reduce_max(x))
    print("Los valores mínimos y máximos después de normalizar son {} y {}.".format(min, max))


# Implementaremos dos tipos diferentes de arquitecturas y las compararemos entre ellas.

num_classes = 1 #la clases son 1 porqué es binario (0 o 1 = effusion o normal)

model = Sequential()

# Add each layer to the model
model.add(Conv2D(128, kernel_size=3,
            activation="relu", input_shape=(64,64,1)))
model.add(MaxPooling2D())
model.add(Conv2D(32, kernel_size=3, activation="relu"))
model.add(MaxPooling2D())
model.add(Conv2D(64, kernel_size=3, activation="relu"))
model.add(MaxPooling2D())
model.add(Flatten())
model.add(Dense(128,activation="relu"))
model.add(Dense(num_classes, activation="sigmoid"))

model.summary()

# Visualizamos el modelo
plot_model(model, show_shapes=True)


# Se cumplen las condiciones requeridas: no hay más de 6 capas convolucionales (3), incluimos pool layers, en la parte inferior las capas totalmente conectadas tienen 128 y 32 nodos, respectivamente, la capa de salida tienen una activación de tipo sigmoide y hay más de 60,000 parametros (351,841).
# Seguidamente compilamos el modelo y lo entrenamos:

model.compile(optimizer='adam',
    loss='binary_crossentropy',
    metrics=['accuracy'])

epochs=30 # entrenamos con 30 epochs

history = model.fit(
                    norm_train,
                    epochs=epochs,
                    validation_data=norm_test,
                    batch_size=batch_size
                )

# Evaluamos el rendimiento usando curvas ROC y las métricas de calidad.

# metricas de calidad
acc = history.history['accuracy']
val_acc = history.history['val_accuracy']

loss = history.history['loss']
val_loss = history.history['val_loss']

epochs_range = range(epochs)

plt.figure(figsize=(8, 8))
plt.subplot(1, 2, 1)
plt.plot(epochs_range, acc, label='Training Accuracy')
plt.plot(epochs_range, val_acc, label='Validation Accuracy')
plt.legend(loc='lower right')
plt.title('Training and Validation Accuracy')

plt.subplot(1, 2, 2)
plt.plot(epochs_range, loss, label='Training Loss')
plt.plot(epochs_range, val_loss, label='Validation Loss')
plt.legend(loc='upper right')
plt.title('Training and Validation Loss')
plt.show()

# curva ROC
# primero realizamos las predicciones de set the testing
y_test_pred = model.predict(norm_test)
# obtenemos valores reales
y_test = y_test = np.concatenate([y for x, y in norm_test], axis=0)


#usamos el paquete sklearn.metrics para calcular los valores de falsos positivos y positivos reales
from sklearn.metrics import roc_curve,roc_auc_score

fpr , tpr , thresholds = roc_curve(y_test , y_test_pred)


def plot_roc_curve(fpr,tpr): 
  plt.plot(fpr,tpr) 
  plt.axis([0,1,0,1]) 
  plt.xlabel('False Positive Rate') 
  plt.ylabel('True Positive Rate')
  plt.title('ROC curve on 1st model') 
  plt.show()    
  
plot_roc_curve (fpr,tpr)

auc_score = roc_auc_score(y_test,y_test_pred)
print("El valor AUC para este modelo es de {}".format(auc_score))


# Vemos que para este modelo la accuracy del training set y del testing set son muy dispares, indicando un sobreentrenamiento, aunque en el caso del testing set se llegue a un buen 0.85. Aún así la curva ROC demuestra un buen valor predictivo del modelo con un valor de AUC de 0.91.
# Asímismo, intentaremos mejorar estos parámetros en el segundo modelo, incluyendo más capas de convuloción y una capa de dropout.

# Creamos, entrenamos y evaluamos el segundo modelo:

# Nuevo modelo
model2 = Sequential()

model2.add(Conv2D(128, kernel_size=3, activation="relu", input_shape=(64,64,1)))
model2.add(MaxPooling2D())
model2.add(Conv2D(32, kernel_size=3, activation="relu"))
model2.add(MaxPooling2D())
model2.add(Conv2D(64, kernel_size=3, activation="relu"))
model2.add(MaxPooling2D())
model2.add(Conv2D(16, kernel_size=3, activation="relu"))
model2.add(MaxPooling2D())
model2.add(Dropout(0.3)) # añadimos una capa de dropout=0.3
model2.add(Flatten())
model2.add(Dense(128,activation="relu"))
model2.add(Dense(num_classes, activation="sigmoid"))

model2.summary()
plot_model(model2, show_shapes=True)

model2.compile(optimizer='adam',
                loss='binary_crossentropy',
                metrics=['accuracy'])

epochs=30 # entrenamos con 30 epochs

history2 = model2.fit(
                    norm_train,
                    epochs=epochs,
                    validation_data=norm_test,
                    batch_size=batch_size
                )

# metricas de calidad
acc = history2.history['accuracy']
val_acc = history2.history['val_accuracy']

loss = history2.history['loss']
val_loss = history2.history['val_loss']

epochs_range = range(epochs)

plt.figure(figsize=(8, 8))
plt.subplot(1, 2, 1)
plt.plot(epochs_range, acc, label='Training Accuracy')
plt.plot(epochs_range, val_acc, label='Validation Accuracy')
plt.legend(loc='lower right')
plt.title('Training and Validation Accuracy')

plt.subplot(1, 2, 2)
plt.plot(epochs_range, loss, label='Training Loss')
plt.plot(epochs_range, val_loss, label='Validation Loss')
plt.legend(loc='upper right')
plt.title('Training and Validation Loss')
plt.show()

# curva ROC
# primero realizamos las predicciones de set the testing
y_test_pred2 = model2.predict(norm_test)
fpr , tpr , thresholds = roc_curve(y_test , y_test_pred2)

def plot_roc_curve(fpr,tpr): 
  plt.plot(fpr,tpr) 
  plt.axis([0,1,0,1]) 
  plt.xlabel('False Positive Rate') 
  plt.ylabel('True Positive Rate')
  plt.title('ROC curve on 2nd model') 
  plt.show()  
plot_roc_curve (fpr,tpr)

auc_score2 = roc_auc_score(y_test,y_test_pred2)
print("El valor AUC para este modelo es de {}".format(auc_score))




