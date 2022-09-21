"""
Name : greeks_calc.py
Author : Meet Kakadiya
Contact : kakadiyameet007.mk@gmail.com
Time    : 05-Jul-21
Desc:
"""

"""
mibian docs
https://github.com/yassinemaaroufi/MibianLib

Documentation
-------------
BS - Black-Scholes        Used for pricing European options on stocks without dividends
BS([underlyingPrice, strikePrice, interestRate, daysToExpiration], volatility=x, callPrice=y, putPrice=z)

eg: 
c = mibian.BS([1.4565, 1.45, 1, 30], volatility=20)
c.callPrice               Returns the call price
c.putPrice                Returns the put price
c.callDelta               Returns the call delta
c.putDelta                Returns the put delta
c.callDelta2              Returns the call dual delta
c.putDelta2               Returns the put dual delta
c.callTheta               Returns the call theta
c.putTheta                Returns the put theta
c.callRho                 Returns the call rho
c.putRho                  Returns the put rho
c.vega                    Returns the option vega
c.gamma                   Returns the option gamma

c = mibian.BS([1.4565, 1.45, 1, 30], callPrice=0.0359)
c.impliedVolatility       Returns the implied volatility from the call price

c = mibian.BS([1.4565, 1.45, 1, 30], putPrice=0.0306)
c.impliedVolatility       Returns the implied volatility from the put price

c = mibian.BS([1.4565, 1.45, 1, 30], callPrice=0.0359, putPrice=0.0306)
c.putCallParity           Returns the put-call parity


GK - Garman-Kohlhagen     Used for pricing European options on currencies
GK([underlyingPrice, strikePrice, domesticRate, foreignRate, daysToExpiration], volatility=x, callPrice=y, putPrice=z)

eg: 
c = mibian.GK([1.4565, 1.45, 1, 2, 30], volatility=20)
c.callPrice               Returns the call price
c.putPrice                Returns the put price
c.callDelta               Returns the call delta
c.putDelta                Returns the put delta
c.callDelta2              Returns the call dual delta
c.putDelta2               Returns the put dual delta
c.callTheta               Returns the call theta
c.putTheta                Returns the put theta
c.callRhoD                Returns the call domestic rho
c.putRhoD                 Returns the put domestic rho
c.callRhoF                Returns the call foreign rho
c.putRhoF                 Returns the call foreign rho
c.vega                    Returns the option vega
c.gamma                   Returns the option gamma

c = mibian.GK([1.4565, 1.45, 1, 2, 30], callPrice=0.0359)
c.impliedVolatility       Returns the implied volatility from the call price

c = mibian.GK([1.4565, 1.45, 1, 2, 30], putPrice=0.03)
c.impliedVolatility       Returns the implied volatility from the put price

c = mibian.GK([1.4565, 1.45, 1, 2, 30], callPrice=0.0359, putPrice=0.03)
c.putCallParity           Returns the put-call parity

Me - Merton               Used for pricing European options on stocks with dividends
Me([underlyingPrice, strikePrice, interestRate, annualDividends, daysToExpiration], volatility=x, callPrice=y, putPrice=z)

eg: 
c = mibian.Me([52, 50, 1, 1, 30], volatility=20)
c.callPrice               Returns the call price
c.putPrice                Returns the put price
c.callDelta               Returns the call delta
c.putDelta                Returns the put delta
c.callDelta2              Returns the call dual delta
c.putDelta2               Returns the put dual delta
c.callTheta               Returns the call theta
c.putTheta                Returns the put theta
c.callRho                 Returns the call rho
c.putRho                  Returns the put rho
c.vega                    Returns the option vega
c.gamma                   Returns the option gamma

c = mibian.Me([52, 50, 1, 1, 30], callPrice=0.0359)
c.impliedVolatility       Returns the implied volatility from the call price

c = mibian.Me([52, 50, 1, 1, 30], putPrice=0.0306)
c.impliedVolatility       Returns the implied volatility from the put price

c = mibian.Me([52, 50, 1, 1, 30], callPrice=0.0359, putPrice=0.0306)
c.putCallParity           Returns the put-call parity



"""


import datetime as dt
import pandas as pd
import numpy as np
from dateutil.relativedelta import relativedelta
import os
import logging
import sys
import pause
import time
import json
from pprint import pprint

import mibian

class Option_greeks:
    """This class gets all the valeue required for calculating options greeks based on paramters"""
    def __init__(self,underlying, strike_price, interest, days_to_expiry
                  ,call_iv=None,put_iv=None,call_price=None ,put_price=None):
        self.underlying=underlying
        self.strike_price= strike_price
        self.interest=interest
        self.days_to_expiry=days_to_expiry
        self.call_iv=call_iv
        self.put_iv=put_iv
        self.call_price=call_price
        self.put_price=put_price


    def get_call_greeks(self):

        obj = mibian.BS([self.underlying, self.strike_price, self.interest, self.days_to_expiry], volatility=  self.call_iv,callPrice=  self.call_price)
        call_greeks = {
            'call_price' : obj.callPrice,
            'call_delta' : obj.callDelta,
            'call_theta' : obj.callTheta,
            'call_iv'    : obj.impliedVolatility,
            'call_rho'   : obj.callRho,
            'call_vega'  : obj.vega,
            'call_gamma' : obj.gamma

        }
        if  self.call_iv!=None:
            call_greeks[ 'call_iv'] = self.call_iv
        if self.call_price !=None:
            call_greeks['call_price'] = self.call_price

        return (call_greeks)

    def get_put_greeks(self):

        obj = mibian.BS([self.underlying, self.strike_price, self.interest, self.days_to_expiry],
                        volatility=self.put_iv, putPrice=self.put_price)
        put_greeks = {
            'put_price': obj.putPrice,
            'put_delta': obj.putDelta,
            'put_theta': obj.putTheta,
            'put_iv': obj.impliedVolatility,
            'put_rho': obj.putRho,
            'put_vega': obj.vega,
            'put_gamma': obj.gamma

        }
        if self.put_iv != None:
            put_greeks['put_iv'] = self.put_iv
        if self.put_price != None:
            put_greeks['put_price'] = self.put_price

        return (put_greeks)


    """
    #RAW BSM 
    
    def black_scholes_merton(stock_price, strike_price, rate, time, volatility, dividend=0.0):

        '''Function that estimates the value of a call and put option using the Black Scholes Merton Model.
      
          Parameters
          ----------
          stock_price: Spot market value of the underlying asset
          strike_price: Strike price of the options contract
          rate: Risk free rate
          time: Time to expiration for the options contract
          volatility: Volatility of the asset
          dividend: Dividend or yield of the asset, with a default value set to zero
        
          Returns
          -------
          [call,put]: Returns a list containing the estimated call and put value of the option contract
          '''
        
          d1 = (log(stock_price/strike_price) + (rate - dividend + volatility**2/2) * time)/(volatility * time**.5)
          d2 = d1 - volatility * time**.5
        
          call = stats.norm.cdf(d1) * stock_price*e**(-dividend*time) - stats.norm.cdf(d2)*strike_price*e**(-rate*time)
          put = stats.norm.cdf(-d2)*strike_price*e**(-rate * time) - stats.norm.cdf(-d1) * stock_price*e**(-dividend*time)
        
          return [call, put]
          
    black_scholes_merton(105,100,.05,1,.25,.01)
    
 
 
 

    """



if __name__=='__main__':
    underlying = 16266.15
    strike_price = 16200
    interest = 10
    days_to_expiry = 6.0
    call_iv = 20.55
    put_iv = 20.55
    call_price =  None
    put_price = None


    option1= Option_greeks(underlying,strike_price ,interest, days_to_expiry,call_iv,put_iv,call_price ,put_price)

    pprint(option1.get_call_greeks())
    print("_____________________________________________________________________________________________________________")
    pprint(option1.get_put_greeks())

