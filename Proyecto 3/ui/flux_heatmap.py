import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
      
#######################################################################
#CLASES
#Clase para hacer zoom en el mapa con scroll
#Obtenido de: https://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
class ZoomPan:
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None


    def zoom_factory(self, ax, base_scale = 2.):
        def zoom(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            xdata = event.xdata 
            ydata = event.ydata 

            if event.button == 'down':
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                scale_factor = base_scale
            else:
                scale_factor = 1
                print event.button

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

            ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            ax.figure.canvas.draw()

        fig = ax.get_figure() 
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, ax):
        def onPress(event):
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            if self.press is None: return
            if event.inaxes != ax: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure() 

        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        return onMotion

#######################################################################
#FUNCIONES
#Lectura y manipulacion de archivos 
def castToNum(pmap):
    for i in range(len(pmap)):
        for j in range(len(pmap[i])):
            if pmap[i][j] != '\n':
                var = float(pmap[i][j])
                pmap[i][j] = var
            else:
                del pmap[i][j]
    return pmap

#Funcion que carga archivo 
def readFile(pfile):
	hmap = []
	try:
		if pfile == "mapa_calor":
			tmp_file = open(pfile + ".txt", "r")
			line = tmp_file.readline()
			row = line.split("#")
			vector = int(row[0]) 
			while line != "":
				hmap += [row]
				line = tmp_file.readline()
				row = line.split("#")
			del hmap[0]
			length = len(hmap)
			hm = castToNum(hmap)
			tmp_file.close()
			show_maps(hm, vector, length)
		else:
			tmp_file = open(pfile + ".txt", "r")
			line = tmp_file.readline()
			row = line.split("#")
			while line != "":
				hmap += [row]
				line = tmp_file.readline()
				row = line.split("#")
			#del hmap[len(hmap)-1]
			tmp_file.close()
			return castToNum(hmap)
	except:
		print ("Hubo algun problema y no se pudo cargar el archivo")
		

def 


#Muestra los mapas de calor y vectorial
def show_maps(R, vector, length):
	R.reverse()
	R = np.array(R)
	color_map = plt.cm.jet
	plt.clf()
	
	#Mapa de calor 
	fig, ax = plt.subplots(num=None, figsize=(18, 12), dpi=80, facecolor='#D4D4D4', edgecolor='k')
	fig.canvas.set_window_title('HEATMAP')
	scale = 1.1
	zoom_obj = ZoomPan()
	fig_zoom = zoom_obj.zoom_factory(ax, base_scale = scale)
	fig_pan = zoom_obj.pan_factory(ax)	
	plt.pcolor(R, cmap = color_map)
	plt.axis([0, length, 0, length]) 
	plt.colorbar().set_label('TEMPERATURE')
	plt.gca().set_axis_bgcolor('black')
	
	#Si se desea ver el mapa vectorial
	if (vector):
		U = readFile("mapa_x")
		V = readFile("mapa_y")
		U = np.array(U)
		V = np.array(V)
		plt.quiver(U, V, alpha=.8)
		plt.quiver(U, V, edgecolor='none', facecolor='none', linewidth=.9)
		plt.xticks(visible = False)
		plt.yticks(visible = False)
	plt.show()

#######################################################################
#EJECUCION
readFile("mapa_calor")
