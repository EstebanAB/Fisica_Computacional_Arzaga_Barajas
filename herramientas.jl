
__precompile__() # Este comando es para que julia precompile el paquete

module herramientas

#Programas método Newton [Tarea 4]:

export metodo_newton
"""documentación del método de Newton"""
function metodo_newton(f,df,x0)
    x=x0
    for i in 1:20
       x=x-f(x)/df(x) 
    end
    return x
end

export metodo_newton_numerico
"""documentación del método de Newton"""
function metodo_newton_numerico(f,x0) #Definimos como al principio
    x = x0
    h = 0.1
    for i in 1:100
        x = x-(f(x))/((f(x+h)-f(x))/(h));     #A partir de aquí, uno se da la idea de que existe un cierto error (Tomar esto como hipotesis)
    end
    return x
end;

export metodo_newton_lista
"""documentación del método de Newton tipo lista"""
function metodo_newton_lista(f,df,inits)
    list=zeros(length(inits))  #Si tenemos n particiones se hace una lista con n elementos
    x = 0.0
    for i in 1:length(inits)   #Aplicamos metodo de newton para cada elemento de nuestro linspace
        x = inits[i]
        for n in 1:200         #Ciclo for, con 200 iteraciones (Este es el epsilon para el ejercicio siguiente)
            x = x-(f(x)/df(x))        
        end
        list[i] = x            #Una lista con iteraciones finales
    end
    return list                #Return hasta que ya haya obtenido las n raices
end;

export metodo_newton_epsilon
"""documentación de metodo de newton epsilon"""
function metodo_newton_epsilon(f,df,inits)
    t = []                                  #Hacemos t un vector vacío
    epsilon = 0.0000001                     #Nuestra epsilon
    list = metodo_newton_lista(f,df,inits)  #Por comodidad, llamamos a la funcion de la lista como lista
    push!(t,list[1])                        #El primer push, para que pueda empezar el ciclo for de i:
    for i in 1:length(t)                    #Por cada elemento de t, hacer esto. Notese que con el push se pueden agregar nuevos elementos
        for n in 1:length(list)             #Ciclo for anidado, ahora sí todos los elementos de la lista de raices
            if abs(t[i]-list[n]) > epsilon  #Comparamos todos los elementos de la lista por solo 1 de t
                push!(t,list[n])            #Si el abs de eso es mayor a epsilon, nos indica que es una raiz diferente.
            end                             #Con el push se agrega un elemento nuevo a t y se repite el primer ciclo for y despues el segundo, de tal manera que empezamos con 1 elemento en t, se compara ese solo elemento con toda la lista: se encuentra una raiz nueva, se agrega a t teniendo 2 elementos y cada uno de esos elementos puede compararse con toda la lista hasta obtener raices diferentes y asi consecutivamente.
        end
        return t                           #Que nos de la lista de raices diferentes                          
    end
end

export metodo_newton_varias_condiciones_iniciales

export integracion_rectangulo
"""documentación de integración por método del rectángulo"""
function integracion_rectangulo(f,a,b,numparticiones)        #Copy-paste de las funciones hechas en tarea 6
    anchorec = (b-a)/numparticiones                     
    s = []                                             
    suma = 0.0                                          
    area = 0.0                                          
    ainit = a                                          
    for i in 1:numparticiones                           
        b = (anchorec)*(i) +ainit                       
        a = ainit+(anchorec)*(i-1)              
        area = (b-a)*(f((b+a)/2))                       
        push!(s,area)                                   
    end
    suma = sum(s,1:numparticiones)                      
    return suma
end;

export integracion_trapecio
"""documentación de integración por método de trapecio"""
function integracion_trapecio(f,a,b,numparticiones)
    anchorec = (b-a)/numparticiones                     
    s = []
    suma = 0.0                    
    area = 0.0                                          
    ainit = a                       
    for i in 1:numparticiones                         
        b = (anchorec)*(i) +ainit                      
        a = ainit+(anchorec)*(i-1)              
        area = (b-a)*((f(a)+f(b))/2)                      
        push!(s,area)
    end
    suma = sum(s,1:numparticiones)
    return suma
end;

export integracion_simpson
"""documentación de integración por método de Simpson"""
function integracion_simpson(f,a,b,numparticiones)
    anchorec = (b-a)/numparticiones                     #practicamente lo mismo que el metodo del rectangulo
    s = []
    suma = 0.0                    
    area = 0.0                                          
    ainit = a                       
    for i in 1:numparticiones                         
        b = (anchorec)*(i) +ainit                      
        a = ainit+(anchorec)*(i-1)              
        area = ((b-a)/6)*((f(a)+4f((a+b)/2)+f(b)))      #p.d. regla de simpson
        push!(s,area)
    end
    suma = sum(s,1:numparticiones)
    return suma
end;

export euler
"""documentación del método de Euler explícito"""
function euler(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x) 
     end
     return listx
    end;

export euler_implicito
"""documentación del método de Euler implicito"""
function euler_implicito(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        t2 = i*h+h
        t_final = h*length(list)
        g(u) = u-x-f(x,t)*h
        x = metodo_newton_numerico(g,x)
        push!(listx,x) 
     end
     return listx
    end;

export euler_pm
"""documentación del método de Euler p.m."""
function euler_pm(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        t2 = i*h+h/2
        x = x + +h*f(x+f(x,t)*h/2,t2)
        push!(listx,x) 
     end
     return listx
    end;

export runge_kutta
"""documentación del método de Runge-Kutta"""
function runge_kutta(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        t2 = i*h+h/2
        #Se actualizan las "kutas":
        k1 = f(x,t)
        k2 = f(x+(h/2)*k1,t2)
        k3 = f(x+(h/2)*k2,t2)
        k4 = f(x+h*k3,t)
        x = x+(h/6)*(k1+2*k2+2*k3+k4) #Por definicion Runge_kutta orden 4
        push!(listx,x) 
     end
     return listx
end;


end
