#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from layer_naive import ProdL, AddL


def apples():
    apple = 100                 # price of an apple
    num_apples  = 2             # number of apples
    tax   = 1.1                 # tax leverage
    prod_apple_layer = ProdL()  
    prod_tax_layer = ProdL()    

    apples_price = prod_apple_layer.forward(apple, num_apples)
    price = prod_tax_layer.forward(apples_price, tax)

    print("money should be payed to shop for {} apples with tax leverage {} and one-apple-price of {} : {}".format(num_apples, tax, apple, price))


    dprice = 1
    dapples_price, dtax = prod_tax_layer.backward(dprice)
    print("derivative of the price w.r.t the price of the apples: {} \nderivative of the price w.r.t the tax leverage: {} yen".format(dapples_price, dtax))
    dapple, dapples_num = prod_apple_layer.backward(dapples_price)
    print("derivative of the price w.r.t the price for each apple: {} \nderivative of the price w.r.t the number of the apples: {} yen".format(dapple, dapples_num))


def apples_and_oranges():
    apple = 100
    orange = 150
    num_apples = 2
    num_oranges = 3
    tax = 1.1

    prod_apple_layer = ProdL()
    prod_orange_layer = ProdL()
    add_apple_orange_layer = AddL()
    prod_tax_layer = ProdL()

    apples_price = prod_apple_layer.forward(apple, num_apples)
    oranges_price = prod_orange_layer.forward(orange, num_oranges)
    apples_oranges_price = add_apple_orange_layer.forward(apples_price, oranges_price)
    price = prod_tax_layer.forward(apples_oranges_price, tax)

    print("money should be paied to shop, calculated by the forward procedures: {}".format(price))

    dprice = 1
    dapples_oranges_price, dtax = prod_tax_layer.backward(dprice)
    print("derivative of the price w.r.t the price of the apples and oranges: {} \nderivative of the price w.r.t the tax leverage: {} yen".format(dapples_oranges_price, dtax))
    dapples_price, doranges_price = add_apple_orange_layer.backward(dapples_oranges_price)
    print("derivative of the price w.r.t the price of the apples: {} \nderivative of the price w.r.t the price of the oranges: {}".format(dapples_price, doranges_price))
    dorange, dnum_oranges = prod_orange_layer.backward(doranges_price)
    print("derivative of the price w.r.t the price for each orange: {} \nderivative of the price w.r.t the number of the oranges: {}".format(dorange, dnum_oranges))
    dapple, dnum_apples = prod_apple_layer.backward(dapples_price)
    print("derivative of the price w.r.t the price for each apple: {} \nderivative of the price w.r.t the number of the apples: {}".format(dapple, dnum_apples))





apples()

apples_and_oranges()
