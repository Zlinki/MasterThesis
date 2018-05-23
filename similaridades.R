setwd("C:\\Users\\Daniel\\ProyectoSenado")
library(data.table)
datos<-fread("reporte_votacioncsv.csv")



#la siguiente funcion permite tratar los dataframes para obtener un dataframe con el que se pueda trabajar


procesarDatos = function(votaciones){
  
  #1 cambia la columna de voto por un factor
  votaciones$voto = as.factor(votaciones$voto)
  #2 cambia la columna de tema por un factor
  votaciones$tema_principal = as.factor(votaciones$tema_principal)
  votaciones$id_votacion= as.factor(votaciones$id_votacion)
  votaciones$voto= as.factor(votaciones$voto)
  
  return(votaciones)
}


#partir los datos para tener senadores y representantes por aparte

datosSenadores = datos[datos$corporacion=="senado",]
datosRepresentantes = datos[datos$corporacion=="camara",]
rm(datos)
#la siguiente funcion calcula el la matriz de likehood entre dos senadores.




#la siguiente funcion calcula, para un tema principal, la similaridad entre las votaciones de dos senadores.
#nota IMPORTANTE: EL DATA frame pasado en votaciones debe ser el dataframe de TODOS  los datos unicamente filtrados
# por senado y camara.


CalcularSimilaridadTema=function(tema,votaciones)
{
  #hace el subset al tema requerido
  votacionesTema = votaciones[votaciones$tema_principal==tema,]
  congresistas =  unique(votaciones$id_congresista)
  ordenados = sort(congresistas)
  similaridad = matrix(0,nrow =length(ordenados),ncol = length(unique(ordenados)))
  totales = matrix(0.00001,nrow =length(ordenados),ncol = length(unique(ordenados)))  
  
  proyectos = unique(votacionesTema$id_proyecto)
for(i in 1:length(ordenados))
 {
   j=1
   for(j in 1:i-1)
     {
        for(k in proyectos)
          {
          proyecto = votacionesTema[votacionesTema$id_proyecto==k,]
          votosi=proyecto[proyecto$id_congresista==i,]
          votosj=proyecto[proyecto$id_congresista==j,]
          
          if(nrow(votosi)>0 && nrow(votosj)>0)
          {
          votoi =names(sort(summary(votosi$voto),decreasing = T)[1])
          votoj =names(sort(summary(votosj$voto),decreasing = T)[1])
          
            if(votoi==votoj)
             {
              similaridad[i,j]=similaridad[i,j]+1
              totales[i,j]=totales[i,j]+1
             }
             if(votoi!=votoj)
             {
               totales[i,j]=totales[i,j]+1
             }
          
            }     
          }
     j=j+1
     }
}
similaridad=similaridad*(1/totales)
similaridad = (similaridad+t(similaridad))  
return(round(similaridad))
#return(list(similaridad,totales))
}


temas = unique(datosSenadores$tema_principal)
primeraRespuesta = CalcularSimilaridadTema("Economia",datosSenadores)
popo = graph_from_adjacency_matrix(primeraRespuesta,mode="undirected")
plot(popo, vertex.size=3, vertex.label=NA)


#funcion para pintar el grafo
pintar= function(tema,votaciones)
{
  resp = CalcularSimilaridadTema(tema,votaciones)
  grafo = graph_from_adjacency_matrix(resp,mode="undirected")
  plot(grafo, vertex.size=3, vertex.label=NA) 
  
  
  
}

