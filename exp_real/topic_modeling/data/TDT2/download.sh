#!/bin/bash

# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/2Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/3Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/4Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/5Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/6Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/7Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/8Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/9Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/10Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/15Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/20Class.zip .
# wget http://www.cad.zju.edu.cn/home/dengcai/Data/TDT2/top30/25Class.zip .

for (( i = 1; i < 11; i++ )); do
    atool -x  "$i"Class.zip
done
