# Physique
## Webscaping (i.e. getting data from the web) reliable websites for Physical constants and making them (easily) accessible locally with (saved) **pandas** dataframes

P.S. I recently discovered juypter notebooks and would point you there, to the [`Physique_jupyter.ipynb`](https://github.com/ernestyalumni/Propulsion/blob/master/Physique/Physique_jupyter.ipynb) jupyter notebook, to see how to use and modify Physique (to webscrape other websites).  Rather in this README.md, I'll point out features that have been added:

* NIST Official Conversions 
* NIST Fundamental Constants = `Physique.FundConst`
* JPL Solar System Dynamics
	- Planets and Pluto

## Using Physique from a working directory not containing Physique itself


```python
import os, sys
```

Get the current directory


```python
currentdir = os.getcwd(); os.getcwd();
```

Then append the directory containing the Physique package/library (it's just a folder) with `sys.path.append`; the absolute path for where I placed it just happened to be `"/Users/ernestyeung/Documents/TeslaModelSP85D"`: substitute that for the absolute path you find (look at your Finder or File organizing program) 


```python
sys.path.append("/Users/ernestyeung/Documents/TeslaModelSP85D")
```


```python
import Physique
```

    Physique
    Physique/rawdata/


Programming note: `__init__.py` in the main directory uses `os.path.dirname(__file__)` with `__file__` (literally that, it's not a placeholder name) being the string with the absolute pathname of the "file from which the module was loaded, if it was loaded from a file" (cf. [stackoverflow  Python __file__ attribute absolute or relative?](http://stackoverflow.com/questions/7116889/python-file-attribute-absolute-or-relative)), i.e. "When a module is loaded in Python, `__file__` is set to its name. You can then use that with other functions to find the directory that the file is located in." (cf. [stackoverflow  what does the __file__ wildcard mean/do?](http://stackoverflow.com/questions/9271464/what-does-the-file-wildcard-mean-do/9271617))  

# NIST Official Conversions

This is the pandas DataFrame containing all the NIST Official Conversions to SI.  


```python
convDF = Physique.conv
```


```python
convDF.columns
```




    Index([u'Toconvertfrom', u'to', u'Multiplyby'], dtype='object')



From the list of columns, search for the quantity you desired by trying out different search terms: e.g. I'm reading Huzel and Huang's **Modern Engineering for Design of Liquid-Propellant Rocket Engines** and I want to know how to convert from
* lb (pound or pound-force) for thrust into force in Newton (N)
* psia (pounds per square inch absolute) for (chamber) pressure into pressure in Pascal (Pa)

We can try to look up the U.S. or Imperial units from the `Toconvertfrom` column.  


```python
convDF[convDF['Toconvertfrom'].str.contains("pound-force ")]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Toconvertfrom</th>
      <th>to</th>
      <th>Multiplyby</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>175</th>
      <td>foot pound-force (ft · lbf)</td>
      <td>joule (J)</td>
      <td>1.355818</td>
    </tr>
    <tr>
      <th>176</th>
      <td>foot pound-force per hour (ft · lbf/h)</td>
      <td>watt (W)</td>
      <td>0.0003766161</td>
    </tr>
    <tr>
      <th>177</th>
      <td>foot pound-force per minute (ft · lbf/min)</td>
      <td>watt (W)</td>
      <td>0.02259697</td>
    </tr>
    <tr>
      <th>178</th>
      <td>foot pound-force per second (ft · lbf/s)</td>
      <td>watt (W)</td>
      <td>1.355818</td>
    </tr>
    <tr>
      <th>340</th>
      <td>pound-force (lbf) 23</td>
      <td>newton (N)</td>
      <td>4.448222</td>
    </tr>
    <tr>
      <th>341</th>
      <td>pound-force foot (lbf · ft)</td>
      <td>newton meter (N · m)</td>
      <td>1.355818</td>
    </tr>
    <tr>
      <th>342</th>
      <td>pound-force foot per inch (lbf · ft/in)</td>
      <td>newton meter per meter (N · m/m)</td>
      <td>53.37866</td>
    </tr>
    <tr>
      <th>343</th>
      <td>pound-force inch (lbf · in)</td>
      <td>newton meter (N · m)</td>
      <td>0.1129848</td>
    </tr>
    <tr>
      <th>344</th>
      <td>pound-force inch per inch (lbf · in/in)</td>
      <td>newton meter per meter (N · m/m)</td>
      <td>4.448222</td>
    </tr>
    <tr>
      <th>345</th>
      <td>pound-force per foot (lbf/ft)</td>
      <td>newton per meter (N/m)</td>
      <td>14.59390</td>
    </tr>
    <tr>
      <th>346</th>
      <td>pound-force per inch (lbf/in)</td>
      <td>newton per meter (N/m)</td>
      <td>175.1268</td>
    </tr>
    <tr>
      <th>347</th>
      <td>pound-force per pound (lbf/lb) \n    (thrust t...</td>
      <td>newton per kilogram (N/kg)</td>
      <td>9.80665</td>
    </tr>
    <tr>
      <th>348</th>
      <td>pound-force per square foot (lbf/ft2)</td>
      <td>pascal (Pa)</td>
      <td>47.88026</td>
    </tr>
    <tr>
      <th>349</th>
      <td>pound-force per square inch (psi) \n    (lbf/in2)</td>
      <td>pascal (Pa)</td>
      <td>6894.757</td>
    </tr>
    <tr>
      <th>350</th>
      <td>pound-force per square inch (psi) (lbf/in2)</td>
      <td>kilopascal (kPa)</td>
      <td>6.894757</td>
    </tr>
    <tr>
      <th>351</th>
      <td>pound-force second per square foot \n    (lbf ...</td>
      <td>pascal second (Pa · s)</td>
      <td>47.88026</td>
    </tr>
    <tr>
      <th>352</th>
      <td>pound-force second per square inch \n    (lbf ...</td>
      <td>pascal second (Pa · s)</td>
      <td>6894.757</td>
    </tr>
    <tr>
      <th>372</th>
      <td>psi (pound-force per square inch) (lbf/in2)</td>
      <td>pascal (Pa)</td>
      <td>6894.757</td>
    </tr>
    <tr>
      <th>373</th>
      <td>psi (pound-force per square inch) (lbf/in2)</td>
      <td>kilopascal (kPa)</td>
      <td>6.894757</td>
    </tr>
  </tbody>
</table>
</div>



Or we can look up the SI unit we want to convert to.


```python
convDF[convDF['to'].str.contains("newton ")]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Toconvertfrom</th>
      <th>to</th>
      <th>Multiplyby</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>137</th>
      <td>dyne (dyn)</td>
      <td>newton (N)</td>
      <td>0.000010</td>
    </tr>
    <tr>
      <th>138</th>
      <td>dyne centimeter (dyn · cm)</td>
      <td>newton meter (N · m)</td>
      <td>1.0E-7</td>
    </tr>
    <tr>
      <th>238</th>
      <td>kilogram-force (kgf)</td>
      <td>newton (N)</td>
      <td>9.80665</td>
    </tr>
    <tr>
      <th>239</th>
      <td>kilogram-force meter (kgf · m)</td>
      <td>newton meter (N · m)</td>
      <td>9.80665</td>
    </tr>
    <tr>
      <th>247</th>
      <td>kilopond (kilogram-force) (kp)</td>
      <td>newton (N)</td>
      <td>9.80665</td>
    </tr>
    <tr>
      <th>250</th>
      <td>kip (1 kip= 1000 lbf)</td>
      <td>newton (N)</td>
      <td>4448.222</td>
    </tr>
    <tr>
      <th>251</th>
      <td>kip (1 kip= 1000 lbf)</td>
      <td>kilonewton (kN)</td>
      <td>4.448222</td>
    </tr>
    <tr>
      <th>300</th>
      <td>ounce (avoirdupois)-force (ozf)</td>
      <td>newton (N)</td>
      <td>0.2780139</td>
    </tr>
    <tr>
      <th>301</th>
      <td>ounce (avoirdupois)-force inch (ozf · in)</td>
      <td>newton meter (N · m)</td>
      <td>0.007061552</td>
    </tr>
    <tr>
      <th>302</th>
      <td>ounce (avoirdupois)-force inch (ozf · in)</td>
      <td>millinewton meter (mN · m)</td>
      <td>7.061552</td>
    </tr>
    <tr>
      <th>336</th>
      <td>poundal</td>
      <td>newton (N)</td>
      <td>0.1382550</td>
    </tr>
    <tr>
      <th>340</th>
      <td>pound-force (lbf) 23</td>
      <td>newton (N)</td>
      <td>4.448222</td>
    </tr>
    <tr>
      <th>341</th>
      <td>pound-force foot (lbf · ft)</td>
      <td>newton meter (N · m)</td>
      <td>1.355818</td>
    </tr>
    <tr>
      <th>342</th>
      <td>pound-force foot per inch (lbf · ft/in)</td>
      <td>newton meter per meter (N · m/m)</td>
      <td>53.37866</td>
    </tr>
    <tr>
      <th>343</th>
      <td>pound-force inch (lbf · in)</td>
      <td>newton meter (N · m)</td>
      <td>0.1129848</td>
    </tr>
    <tr>
      <th>344</th>
      <td>pound-force inch per inch (lbf · in/in)</td>
      <td>newton meter per meter (N · m/m)</td>
      <td>4.448222</td>
    </tr>
    <tr>
      <th>345</th>
      <td>pound-force per foot (lbf/ft)</td>
      <td>newton per meter (N/m)</td>
      <td>14.59390</td>
    </tr>
    <tr>
      <th>346</th>
      <td>pound-force per inch (lbf/in)</td>
      <td>newton per meter (N/m)</td>
      <td>175.1268</td>
    </tr>
    <tr>
      <th>347</th>
      <td>pound-force per pound (lbf/lb) \n    (thrust t...</td>
      <td>newton per kilogram (N/kg)</td>
      <td>9.80665</td>
    </tr>
    <tr>
      <th>423</th>
      <td>ton-force (2000 lbf)</td>
      <td>newton (N)</td>
      <td>8896.443</td>
    </tr>
    <tr>
      <th>424</th>
      <td>ton-force (2000 lbf)</td>
      <td>kilonewton (kN)</td>
      <td>8.896443</td>
    </tr>
  </tbody>
</table>
</div>



Look at what you want and see the index; it happens to be 340 in this example.  


```python
lbf2N = convDF.loc[340,:]; lbf2N
```




    Toconvertfrom    pound-force (lbf) 23
    to                         newton (N)
    Multiplyby                   4.448222
    Name: 340, dtype: object



Then the attributes can accessed by the column names.  


```python
print lbf2N.Toconvertfrom, lbf2N.to, lbf2N.Multiplyby
```

    pound-force (lbf) 23 newton (N) 4.448222


So for example, the reusable SSME delivers a vacuum thrust of 470000 lb or 


```python
print 470000*lbf2N.Multiplyby, lbf2N.to 
```

    2090664.340000 newton (N)


To obtain the conversion for pressure in psia, which we search for with "psi" 


```python
convDF[convDF['Toconvertfrom'].str.match("psi")]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Toconvertfrom</th>
      <th>to</th>
      <th>Multiplyby</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>372</th>
      <td>psi (pound-force per square inch) (lbf/in2)</td>
      <td>pascal (Pa)</td>
      <td>6894.757</td>
    </tr>
    <tr>
      <th>373</th>
      <td>psi (pound-force per square inch) (lbf/in2)</td>
      <td>kilopascal (kPa)</td>
      <td>6.894757</td>
    </tr>
  </tbody>
</table>
</div>



So for a chamber pressure of 3028 psia for the SSME, 


```python
psi2Pa = convDF.loc[372,:]
```


```python
print 3028*psi2Pa.Multiplyby, psi2Pa.to
```

    20877324.196 pascal (Pa)


Also, get the conversion for atmospheres (atm):


```python
convDF[convDF['Toconvertfrom'].str.match("atm")]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Toconvertfrom</th>
      <th>to</th>
      <th>Multiplyby</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>15</th>
      <td>atmosphere, standard (atm)</td>
      <td>pascal (Pa)</td>
      <td>101325</td>
    </tr>
    <tr>
      <th>16</th>
      <td>atmosphere, standard (atm)</td>
      <td>kilopascal (kPa)</td>
      <td>101.325</td>
    </tr>
    <tr>
      <th>17</th>
      <td>atmosphere, technical (at) 8</td>
      <td>pascal (Pa)</td>
      <td>98066.5</td>
    </tr>
    <tr>
      <th>18</th>
      <td>atmosphere, technical (at) 8</td>
      <td>kilopascal (kPa)</td>
      <td>98.0665</td>
    </tr>
  </tbody>
</table>
</div>




```python
atm2Pa = convDF.loc[15,:]
```


```python
print 3028*psi2Pa.Multiplyby/atm2Pa.Multiplyby, atm2Pa.Toconvertfrom
```

    206.0431699580557611645694547 atmosphere, standard (atm)



```python

```
