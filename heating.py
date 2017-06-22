# -*- coding: utf-8 -*-
#######
# Flight Path Heating
# Author : Susumu Tanaka
# Ref.
# * 宇宙飛行体の熱気体力学
# * Heat Transfer to Satellite Vehicles Re-entering the Atomosphere
# * 超軌道足度飛行体の輻射加熱環境に関する研究

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import environment as env

class heating:
  def __init__(self):
    ######## input ############
    self.logfile_name = 'Momo_dynamics_PlanB_20170619(Rev6.2).csv'
    self.log_all = np.loadtxt(self.logfile_name, delimiter=',', skiprows=1)
    self.time = self.log_all[:,0]
    self.altitude = self.log_all[:,5]
    self.mach = self.log_all[:,22]
    self.array_length = len(self.time)

    self.T_surface_init = 15.0 # [degC] flight object initial temperature
    self.T_surface_init += 273.15 # [K]
    self.R_nose = 0.02 # [m] nose-cone blunt radius
    self.rho_nose = 1270.0 # [kg/m^3] nose-cone material dencity
    self.thickness = 0.025 # [m] nose-cone thickness at stagnation point
    self.specific_heat = 1591.0 # [-] nose-cone material
    self.epsilon = 0.8 # [-] 表面放射率
    ###########################

    self.Re = 6371000 #[m] earth surface
    self.Cp = 1006.0 # [J/kg-K] specific heat at pressure constant of air
    self.sigma = 5.669 * 10**(-8) # Stefan-Boltzmann constant

    self.T0, self.rho0, self.Cs0 = env.std_atmo(0.0) # [K, kg/m^3]
    self.g0 = env.gravity(0.0) # [m/s^2]       


  def calculation(self):
    self.T = np.empty(self.array_length)
    self.rho = np.empty(self.array_length)
    self.Cs = np.empty(self.array_length)
    self.g = np.empty(self.array_length)
    self.vel = np.empty(self.array_length)
    self.R = np.empty(self.array_length)
    self.uc = np.empty(self.array_length)
    self.qconv = np.empty(self.array_length) # convection heating
    self.qrad = np.empty(self.array_length) # radiative heating
    self.T_surface_thinskin = np.empty(self.array_length) # thin-skin method

    def q_convection(R_nose, rho, rho0, vel, uc):
      q_convection = 11030.0 / np.sqrt(R_nose) * (rho / rho0)**0.5 * (np.abs(vel) / uc)**3.05 * 10**4 # [W/m^2]
      return q_convection
    
    def q_radiation(R_nose, vel, rho):
      # ref. Tauberの経験式
      def radiative_heating_velocity_function(vel):
        # input:[m/s]
        vel = np.abs(vel) / 1000.0 # [km/s]
        vel_array = np.array([9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0])
        return_array = np.array([1.5, 4.3, 9.7, 18.5, 35.0, 55.0, 81.0, 115.0, 151.0, 238.0, 359.0, 495.0, 660.0, 850.0, 1065.0, 1313.0, 1550.0, 1780.0, 2040.0])
        inter = interpolate.interp1d(vel_array, return_array, bounds_error = False, fill_value = (return_array[0], return_array[-1]))
        return inter(vel)
      def coef_a(R_nose, vel, rho):
        # input:[m, m/s, kg/m^3]
        a = 1.072 * 10.0**6 * np.abs(vel)**(-1.88) * rho**(-0.325)
        if R_nose <= 1.0:
          return a
        elif R_nose >= 2.0:
          return min(0.5, a)
        else:
          return min(0.6, a)
      q_radiation = 4.736 * 10**4 * R_nose**coef_a(R_nose, vel, rho) * rho**1.22 * radiative_heating_velocity_function(vel) * 10**4 # [W/m^2]
      return q_radiation
    
    def deltaT_surface_thinskin(T_surface_pre, qconv, qrad, sigma, epsilon, specific_heat, rho_nose, thickness):
      return (qconv + qrad - sigma * epsilon * T_surface_pre**4) / (specific_heat * rho_nose * thickness) # [K]
    
    self.T[0], self.rho[0], self.Cs[0] = env.std_atmo(self.altitude[0])
    self.g[0] = env.gravity(self.altitude[0])
    self.vel[0] = self.mach[0] * self.Cs[0]
    self.R[0] = self.Re + self.altitude[0] # [m]
    self.uc[0] = np.sqrt(self.g0 * self.Re**2 / self.R[0]) # [m/s]
    self.qconv[0] = q_convection(self.R_nose, self.rho[0], self.rho0, self.vel[0], self.uc[0]) # convection heating
    self.qrad[0] = q_radiation(self.R_nose, self.vel[0], self.rho[0]) # radiative heating
    self.T_surface_thinskin[0] = self.T_surface_init

    for i in range(1, self.array_length):
      self.dt = self.time[i] - self.time[i-1]
      self.T[i], self.rho[i], self.Cs[i] = env.std_atmo(self.altitude[i])
      self.g[i] = env.gravity(self.altitude[i])
      self.vel[i] = self.mach[i] * self.Cs[i]
      self.R[i] = self.Re + self.altitude[i] # [m]
      self.uc[i] = np.sqrt(self.g0 * self.Re**2 / self.R[i]) # [m/s]
      self.qconv[i] = q_convection(self.R_nose, self.rho[i], self.rho0, self.vel[i], self.uc[i]) # convection heating
      self.qrad[i] = q_radiation(self.R_nose, self.vel[i], self.rho[i]) # radiative heating
      self.T_surface_thinskin[i] = self.T_surface_thinskin[i-1] + self.dt * deltaT_surface_thinskin(self.T_surface_thinskin[i-1], self.qconv[i], self.qrad[i], self.sigma, self.epsilon, self.specific_heat, self.rho_nose, self.thickness) # thin-skin method
      
      print(self.time[i], self.T[i],self.qconv[i], self.qrad[i], self.T_surface_thinskin[i])

if __name__ == '__main__':
  obj = heating()
  obj.calculation()

  # result = np.c_[obj.K, obj.hsw, obj.qconv, obj.qrad, obj.T_surface_radequ, obj.T_surface_thinskin]
  # np.savetxt('log.csv',result, delimiter=',')

  plt.figure(0)
  plt.plot(obj.time, obj.qconv/10**6)
  plt.plot(obj.time, obj.qrad/10**6)
  plt.xlabel('time [sec]')
  plt.ylabel('q_dot [MW/m2]')
  plt.figure(1)
  plt.plot(obj.time, obj.T_surface_thinskin)
  plt.xlabel('time [sec]')
  plt.ylabel('T_surface [K]')

  plt.figure(4)
  plt.plot(obj.time, obj.altitude)
  plt.xlabel('time [sec]')
  plt.ylabel('altitude [m]')
  plt.figure(5)
  plt.plot(obj.time, obj.vel)
  plt.xlabel('time [sec]')
  plt.ylabel('velocity [m/s]')

  plt.show()