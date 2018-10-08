from __future__ import print_function
import sys
import numpy as np
from numpy import *
import tkinter as tk
import math as m 
from visual import *

#Function that normalizes u
def norm(u):
	if u >= 0:
		return 2*pi * (u/(2*pi) - int(u/(2*pi)))
	else:
		return 2*pi + ( 2*pi * (u/(2*pi) - int(u/(2*pi))) )

#Function that calculates the fi function for the N-R method
def fi(epsilon, xi, u):
	return (epsilon*(sin(u)-u*cos(u))+xi)/(1-epsilon*cos(u))

#Function that implements the N-R method for u calculations
def NR (epsilon, xi):
	u_ant = pi
	u_act = fi(epsilon, xi, u_ant)

	while abs(u_act - u_ant) > prec:
		u_ant = u_act
		u_act = fi(epsilon, xi, u_ant)

	return float(norm(u_act))

class Planet:
	def __init__(self, name, a, epsilon, p, i, Omega, bar_omega, static_ball, orbit_drawer, label, t0 = 0):
		self.name = name
		self.a = a
		self.b = a*sqrt(1-epsilon*epsilon)
		self.epsilon = epsilon
		self.p = p
		self.sqrtmu = 2*pi*pow(self.a, 1.5)/self.p #The sqrt of the mu parameter
		self.t0 = t0 #Days from 1/1/1900 until first step at perihelion
		self.static_ball = static_ball
		self.orbit_drawer = orbit_drawer
		self.label = label
		self.i = i
		self.Omega = Omega
		self.bar_omega = bar_omega
		self.omega = bar_omega - Omega

	#Function that calculates xi for the Kepler's equation for a given time
	def calc_xi (self, time):
		return 2*pi / self.p * (time - self.t0)

	#Function that calculates the eccentric anomaly via Newton-Raphson method for a given time
	def u (self, time):
		return NR(self.epsilon, self.calc_xi(time))

	#Function that calculates the derivative of the eccentric anomaly for a given time
	def devu (self, time):
		u = self.u(time)
		return self.sqrtmu / (pow(self.a, 1.5)*(1-self.epsilon*cos(u)))

	#Function that calculates the planet's position for a given time
	def getPos (self, time):
		u = self.u(time)
		row1 = (float(cos(self.Omega)*cos(self.i)*cos(self.omega)+sin(self.Omega)*sin(self.omega)), float(-cos(self.i)*sin(self.omega)*cos(self.Omega)-sin(self.Omega)*cos(self.omega)), float(cos(self.Omega)*sin(self.i)))
		row2 = (float(-cos(self.i)*cos(self.omega)*sin(self.Omega)-cos(self.Omega)*sin(self.omega)), float(sin(self.Omega)*cos(self.i)*sin(self.omega)+cos(self.Omega)*cos(self.omega)), float(sin(self.Omega)*sin(self.i)))
		row3 = (float(sin(self.i)*cos(self.omega)),float(sin(self.i)*sin(self.omega)),float(cos(self.i)))

		basic_pos = ( float(self.a*(cos(u)-self.epsilon)),float(self.a*sqrt(1-self.epsilon*self.epsilon)*sin(u)), 0)

		pos = (float(np.dot(row1, basic_pos)), float(np.dot(row2, basic_pos)), float(np.dot(row3, basic_pos)))
		return pos

	#Function that modifies the planet's ball position
	def setPos (self, time):
		self.static_ball.pos = self.getPos(time)
		self.label.pos = self.static_ball.pos

	#Function that calculates the angular momentum for a given time
	def getAngMomentum (self, time):
		u = self.u(time)
		du = self.devu(time)
		return (0, 0, float(self.a*self.a*sqrt(1-self.epsilon*self.epsilon)*du*(1-self.epsilon*cos(u))) )

	#Function that calculates the energy for a given time
	def getEnergy (self, time):
		u = self.u(time)
		du = self.devu(time)
		return 0.5*self.a*self.a*du*du*(1-self.epsilon*self.epsilon*cos(u)*cos(u)) - self.sqrtmu*self.sqrtmu/(self.a*(1-self.epsilon*cos(u)))

	#Function that calculates the t for a given eccentric anomaly
	def getTimeFromu (self, u):
		return self.p / (2*pi) * (u - self.epsilon*sin(u)) + self.t0

	#Function that draws the planet's orbit
	def draw_orbit(self):
		self.orbit_drawer.trail_object.visible = True
		times = np.linspace(0, self.p, 2000, endpoint=True)
		for t in times:
			self.orbit_drawer.pos = self.getPos(t)

	#Function that draws the planet's orbit and its position for a given time
	def draw(self, time):
		self.static_ball.visible = True
		self.label.visible = True
		self.static_ball.pos = self.getPos(time)
		self.label.pos = self.static_ball.pos
		self.draw_orbit()

	#Function that hide the planet
	def hide(self):
		self.static_ball.visible = False
		self.label.visible = False
		self.orbit_drawer.trail_object.visible = False

	#Function that prints several planets' data for a given time
	def printData (self, time, u):
		print(self.name)
		print("a = ", self.a, " b = ", self.b)
		print("Para t =", time, "dias terrestres.")
		print("Posicion:", self.getPos(time))
		print("Momento angular:", self.getAngMomentum(time))
		print("Anomalia excentrica con N-R:", self.u(time))
		print("Energia:", self.getEnergy(time))
		print("Para el u introducido el t era:", self.getTimeFromu(u))

def planet():
	global t_animation
	global first_time
	t = float(var_text_time.get())
	u = float(var_text_u.get())

	t_animation = t

	if first_time:
		first_time = False
		animate()

	#For the selected planets, draw its position and print several data
	for i in range(0,8):
		if (planet_list[i].get() == 1):
			planets[i].draw(t)
			planets[i].printData(t, u)
		else:
			planets[i].hide();

def playPause():
	global run_animation
	run_animation = not run_animation

def animate():
	global t_animation

	if run_animation:
		t_animation = t_animation + 0.1

		for i in range(0,8):
			if (planet_list[i].get() == 1):
				planets[i].setPos(t_animation)

	master.after(2, animate)


def quit():
	exit()
#Create GUI
master = tk.Tk()

tk.Label(master, text="Selecciona los planetas que desees:").grid(row=0, sticky=tk.W)
var1 = tk.IntVar()
tk.Checkbutton(master, text="Mercurio", variable=var1).grid(row=1, sticky=tk.W)
var2 = tk.IntVar()
tk.Checkbutton(master, text="Venus", variable=var2).grid(row=2, sticky=tk.W)
var3 = tk.IntVar()
tk.Checkbutton(master, text="Tierra", variable=var3).grid(row=3, sticky=tk.W)
var4 = tk.IntVar()
tk.Checkbutton(master, text="Marte", variable=var4).grid(row=4, sticky=tk.W)
var5 = tk.IntVar()
tk.Checkbutton(master, text="Jupiter", variable=var5).grid(row=5, sticky=tk.W)
var6 = tk.IntVar()
tk.Checkbutton(master, text="Saturno", variable=var6).grid(row=6, sticky=tk.W)
var7 = tk.IntVar()
tk.Checkbutton(master, text="Urano", variable=var7).grid(row=7, sticky=tk.W)
var8 = tk.IntVar()
tk.Checkbutton(master, text="Neptuno", variable=var8).grid(row=8, sticky=tk.W)

tk.Label(master, text="Introduce el t (en dias terrestres a partir del 1-enero-1900) para el que quieres obtener la posicion de los planetas: ").grid(row=9)
var_text_time = tk.StringVar()
text_box_time = tk.Entry(master, textvariable = var_text_time).grid(row=9, column=1)

tk.Label(master, text="Introduce el u para el que quieres obtener el tiempo en el que se dio: ").grid(row=10)
var_text_u = tk.StringVar()
text_box_u = tk.Entry(master, textvariable = var_text_u).grid(row=10, column=1)

tk.Button(master, text='Show', command = planet).grid(row=11, sticky=tk.W, pady=4)
tk.Button(master, text='Play/Pause', command = playPause).grid(row=12, sticky=tk.W, pady=4)
tk.Button(master, text='Quit', command = quit).grid(row=13, sticky=tk.W, pady=4)


#Newton-Raphson accuracy
prec = 0.0001

#Global variables for animation
run_animation = False
t_animation = 0
first_time = True

#Create scene and Sun:
scene = display(title='Sistema Solar', x=0, y=0, width=10000, height=10000)
Sun = sphere(pos=(0,0,0), radius=0.03500, color=color.yellow)

#Create planets' balls:
Mercury_static_ball = sphere(radius=0.015,material=materials.wood, visible = False)
Venus_static_ball = sphere(radius=0.015,material=materials.rough, color=color.yellow, visible = False)
Earth_static_ball = sphere(radius=0.0190, material=materials.earth, visible = False)
Mars_static_ball = sphere(radius=0.0190, material=materials.shiny,color=color.orange, visible = False)
Jupiter_static_ball = sphere(radius=0.2, material=materials.marble, visible = False)
Saturn_static_ball = sphere(radius=0.2, material=materials.marble,color=color.orange, visible = False)
Uranus_static_ball = sphere(radius=0.2, material=materials.glass,color=color.green, visible = False)
Neptune_static_ball = sphere(radius=0.2,material=materials.glass, color=color.blue, visible = False)

#Create orbit drawers:
Mercury_orbit_drawer = sphere(radius=0.01,material=materials.wood, make_trail =true, visible = False)
Venus_orbit_drawer = sphere(radius=0.01,material=materials.rough, color=color.yellow, make_trail =true, visible = False)
Earth_orbit_drawer = sphere(radius=0.01, material=materials.earth,make_trail=true, visible = False)
Mars_orbit_drawer = sphere(radius=0.01, material=materials.shiny,color=color.orange, make_trail =true, visible = False)
Jupiter_orbit_drawer = sphere(radius=0.01, material=materials.marble, make_trail =true, visible = False)
Saturn_orbit_drawer = sphere(radius=0.01, material=materials.marble,color=color.orange, make_trail =true, visible = False)
Uranus_orbit_drawer = sphere(radius=0.01, material=materials.glass,color=color.green, make_trail =true, visible = False)
Neptune_orbit_drawer = sphere(radius=0.01,material=materials.glass, color=color.blue, make_trail =true, visible = False)

#Create planets' labels:

Neptune_label = label(text = "Neptuno", color= color.white, visible = False)
Uranus_label = label(text = "Urano", color= color.white, visible = False)
Saturn_label = label(text = "Saturno", color= color.white, visible = False)
Jupiter_label = label(text = "Jupiter", color= color.white, visible = False)
Mars_label = label(text = "Marte", color= color.white, visible = False)
Earth_label = label(text = "Tierra", color= color.white, visible = False)
Venus_label = label(text = "Venus", color= color.white, visible = False)
Mercury_label = label(text = "Mercurio", color= color.white, visible = False)

#Create the planets:
grad2radConst = pi / 180
Mercury = Planet("Mercurio", 0.387, 0.206, 87.97,  grad2radConst * 7, grad2radConst * 47.14, grad2radConst *75.9, Mercury_static_ball, Mercury_orbit_drawer, Mercury_label, 61)
Venus = Planet("Venus", 0.723, 0.007, 224.7, grad2radConst * 3.59, grad2radConst * 75.78, grad2radConst * 130.15, Venus_static_ball, Venus_orbit_drawer, Venus_label, 90)
Earth = Planet("Tierra", 1.0, 0.017, 365.26, 0, 0, grad2radConst * 101.22, Earth_static_ball, Earth_orbit_drawer, Earth_label)
Mars = Planet("Marte", 1.524, 0.093, 686.98, grad2radConst * 1.85, grad2radConst * 48.78, grad2radConst * 334.22, Mars_static_ball, Mars_orbit_drawer, Mars_label, 76)
Jupiter = Planet("Jupiter", 5.203, 0.048, 4332.6, grad2radConst * 1.31, grad2radConst * 99.44, grad2radConst * 12.72, Jupiter_static_ball, Jupiter_orbit_drawer, Jupiter_label, 1612)
Saturn = Planet("Saturno", 9.546, 0.056, 10759, grad2radConst * 2.5, grad2radConst * 112.79, grad2radConst * 91.09, Saturn_static_ball, Saturn_orbit_drawer, Saturn_label, 5528)
Uranus = Planet("Urano", 19.20, 0.047, 30687, grad2radConst * 0.77, grad2radConst * 73.48, grad2radConst * 169.05, Uranus_static_ball, Uranus_orbit_drawer, Uranus_label, 24245)
Neptune = Planet("Neptuno", 30.09, 0.009, 60784, grad2radConst * 1.78, grad2radConst * 130.68, grad2radConst * 43.83, Neptune_static_ball, Neptune_orbit_drawer, Neptune_label, 52122)

planets = []
planets.append(Mercury)
planets.append(Venus)
planets.append(Earth)
planets.append(Mars)
planets.append(Jupiter)
planets.append(Saturn)
planets.append(Uranus)
planets.append(Neptune)


planet_list = []
planet_list.append(var1)
planet_list.append(var2)
planet_list.append(var3)
planet_list.append(var4)
planet_list.append(var5)
planet_list.append(var6)
planet_list.append(var7)
planet_list.append(var8)

#Start Tkinter loop
tk.mainloop()
