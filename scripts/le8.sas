/*LE8*/
/*第一步 将变量中的coding进行重新赋值 555=0.5 444=0.25 300=3*/
data diet24op;
set wmy.diet24allnew;
array alln _numeric_;
do over alln;
if  alln=555 then alln=0.5;
if  alln=444 then alln=0.25;
if  alln=200 then alln=2;
if  alln=300 then alln=3;
if  alln=400 then alln=4;
if  alln=500 then alln=5;
if  alln=600 then alln=6;
end;
run;
/*第二步 处理atypical diet*/
data diet24opty;
set diet24op;
array p(*)  n_104000_0_0 n_104010_0_0 n_104020_0_0 n_104030_0_0 n_104040_0_0 n_104050_0_0 n_104060_0_0 n_104070_0_0 n_104080_0_0 n_104090_0_0 
n_104100_0_0 n_104110_0_0 n_104120_0_0 n_104130_0_0 n_104140_0_0 n_104150_0_0 n_104160_0_0 n_104170_0_0 n_104180_0_0 n_104190_0_0 n_104200_0_0 n_104210_0_0 
n_104220_0_0 n_104230_0_0 n_104240_0_0 n_104250_0_0 n_104260_0_0 n_104270_0_0 n_104280_0_0 n_104290_0_0 n_104300_0_0 n_104310_0_0 n_104320_0_0 n_104330_0_0 
n_104340_0_0 n_104350_0_0 n_104360_0_0 n_104370_0_0 n_104380_0_0 n_104400_0_0 n_104410_0_0 n_104420_0_0 n_104430_0_0 n_104440_0_0 n_104450_0_0 n_104460_0_0 
n_104470_0_0 n_104480_0_0 n_104490_0_0 n_104500_0_0 n_104510_0_0 n_104520_0_0 n_104530_0_0 n_104540_0_0 n_104550_0_0 n_104560_0_0 n_104570_0_0 n_104580_0_0 
n_104590_0_0 n_102410_0_0 n_102420_0_0 n_102430_0_0 n_102440_0_0 n_102450_0_0 n_103270_0_0 n_104000_0_0 n_104010_0_0 n_104110_0_0 n_104120_0_0 n_104280_0_0 
n_100770_0_0 n_100800_0_0 n_100810_0_0 n_100820_0_0 n_100830_0_0 n_100840_0_0 n_100850_0_0 n_100860_0_0 n_102740_0_0 n_102780_0_0 n_102770_0_0 
n_102080_0_0 n_102090_0_0 n_20106_0_0  n_102810_0_0 n_102850_0_0 n_102870_0_0 n_100520_0_0 n_100920_0_0 n_100160_0_0 n_100170_0_0 n_100180_0_0 
n_100190_0_0 n_100200_0_0 n_100210_0_0 n_100220_0_0 n_100230_0_0 n_100370_0_0 n_100380_0_0 n_100490_0_0 n_100500_0_0 n_100540_0_0 n_100550_0_0 n_103000_0_0 
n_103010_0_0 n_103020_0_0 n_103030_0_0 n_103040_0_0 n_103070_0_0 n_103080_0_0 n_103090_0_0 n_30530_0_0  
n_103990_0_0 n_104400_0_0 n_102400_0_0 n_103310_0_0	n_20088_0_0	 n_103280_0_0 n_103260_0_0 n_103290_0_0 n_100760_0_0
n_100940_0_0 n_20091_0_0  n_100950_0_0 n_20092_0_0	n_101020_0_0 n_20093_0_0 n_101090_0_0 n_20094_0_0 n_101160_0_0 
n_102720_0_0 n_101250_0_0 n_101260_0_0 n_101270_0_0 n_102800_0_0 n_102970_0_0; 
do i=1 to dim(p);
if n_20085_0_0 = 3 or   n_20085_0_0 = 4 or  n_20085_0_0 = 6  then p(i)=.;
end;
array f(*)  n_104000_1_0 n_104010_1_0 n_104020_1_0 n_104030_1_0 n_104040_1_0 n_104050_1_0 n_104060_1_0 n_104070_1_0 n_104080_1_0 n_104090_1_0 
n_104100_1_0 n_104110_1_0 n_104120_1_0 n_104130_1_0 n_104140_1_0 n_104150_1_0 n_104160_1_0 n_104170_1_0 n_104180_1_0 n_104190_1_0 n_104200_1_0 n_104210_1_0 
n_104220_1_0 n_104230_1_0 n_104240_1_0 n_104250_1_0 n_104260_1_0 n_104270_1_0 n_104280_1_0 n_104290_1_0 n_104300_1_0 n_104310_1_0 n_104320_1_0 n_104330_1_0 
n_104340_1_0 n_104350_1_0 n_104360_1_0 n_104370_1_0 n_104380_1_0 n_104400_1_0 n_104410_1_0 n_104420_1_0 n_104430_1_0 n_104440_1_0 n_104450_1_0 n_104460_1_0 
n_104470_1_0 n_104480_1_0 n_104490_1_0 n_104500_1_0 n_104510_1_0 n_104520_1_0 n_104530_1_0 n_104540_1_0 n_104550_1_0 n_104560_1_0 n_104570_1_0 n_104580_1_0 
n_104590_1_0 n_102410_1_0 n_102420_1_0 n_102430_1_0 n_102440_1_0 n_102450_1_0 n_103270_1_0 n_104000_1_0 n_104010_1_0 n_104110_1_0 n_104120_1_0 n_104280_1_0 
n_100770_1_0 n_100800_1_0 n_100810_1_0 n_100820_1_0 n_100830_1_0 n_100840_1_0 n_100850_1_0 n_100860_1_0 n_102740_1_0 n_102780_1_0 n_102770_1_0 
n_102080_1_0 n_102090_1_0 n_20106_1_0  n_102810_1_0 n_102850_1_0 n_102870_1_0 n_100520_1_0 n_100920_1_0 n_100160_1_0 n_100170_1_0 n_100180_1_0 
n_100190_1_0 n_100200_1_0 n_100210_1_0 n_100220_1_0 n_100230_1_0 n_100370_1_0 n_100380_1_0 n_100490_1_0 n_100500_1_0 n_100540_1_0 n_100550_1_0 n_103000_1_0 
n_103010_1_0 n_103020_1_0 n_103030_1_0 n_103040_1_0 n_103070_1_0 n_103080_1_0 n_103090_1_0 n_30530_1_0 
n_103990_1_0 n_104400_1_0 n_102400_1_0 n_103310_1_0	n_20088_1_0	 n_103280_1_0 n_103260_1_0 n_103290_1_0 n_100760_1_0
n_100940_1_0 n_20091_1_0  n_100950_1_0 n_20092_1_0	n_101020_1_0 n_20093_1_0 n_101090_1_0 n_20094_1_0 n_101160_1_0 
n_102720_1_0 n_101250_1_0 n_101260_1_0 n_101270_1_0 n_102800_1_0 n_102970_1_0; 
do i=1 to dim(f);
if n_20085_1_0 = 3 or   n_20085_1_0 = 4 or  n_20085_1_0 = 6  then f(i)=.;
end;
array s(*)  n_104000_2_0 n_104010_2_0 n_104020_2_0 n_104030_2_0 n_104040_2_0 n_104050_2_0 n_104060_2_0 n_104070_2_0 n_104080_2_0 n_104090_2_0 
n_104100_2_0 n_104110_2_0 n_104120_2_0 n_104130_2_0 n_104140_2_0 n_104150_2_0 n_104160_2_0 n_104170_2_0 n_104180_2_0 n_104190_2_0 n_104200_2_0 n_104210_2_0 
n_104220_2_0 n_104230_2_0 n_104240_2_0 n_104250_2_0 n_104260_2_0 n_104270_2_0 n_104280_2_0 n_104290_2_0 n_104300_2_0 n_104310_2_0 n_104320_2_0 n_104330_2_0 
n_104340_2_0 n_104350_2_0 n_104360_2_0 n_104370_2_0 n_104380_2_0 n_104400_2_0 n_104410_2_0 n_104420_2_0 n_104430_2_0 n_104440_2_0 n_104450_2_0 n_104460_2_0 
n_104470_2_0 n_104480_2_0 n_104490_2_0 n_104500_2_0 n_104510_2_0 n_104520_2_0 n_104530_2_0 n_104540_2_0 n_104550_2_0 n_104560_2_0 n_104570_2_0 n_104580_2_0 
n_104590_2_0 n_102410_2_0 n_102420_2_0 n_102430_2_0 n_102440_2_0 n_102450_2_0 n_103270_2_0 n_104000_2_0 n_104010_2_0 n_104110_2_0 n_104120_2_0 n_104280_2_0 
 n_100770_2_0 n_100800_2_0 n_100810_2_0 n_100820_2_0 n_100830_2_0 n_100840_2_0 n_100850_2_0 n_100860_2_0 n_102740_2_0 n_102780_2_0 n_102770_2_0 
n_102080_2_0 n_102090_2_0 n_20106_2_0  n_102810_2_0 n_102850_2_0 n_102870_2_0 n_100520_2_0 n_100920_2_0 n_100160_2_0 n_100170_2_0 n_100180_2_0 
n_100190_2_0 n_100200_2_0 n_100210_2_0 n_100220_2_0 n_100230_2_0 n_100370_2_0 n_100380_2_0 n_100490_2_0 n_100500_2_0 n_100540_2_0 n_100550_2_0 n_103000_2_0 
n_103010_2_0 n_103020_2_0 n_103030_2_0 n_103040_2_0 n_103070_2_0 n_103080_2_0 n_103090_2_0 n_30530_2_0 
n_103990_2_0 n_104400_2_0 n_102400_2_0 n_103310_2_0	n_20088_2_0	 n_103280_2_0 n_103260_2_0 n_103290_2_0 n_100760_2_0
n_100940_2_0 n_20091_2_0  n_100950_2_0 n_20092_2_0	n_101020_2_0 n_20093_2_0 n_101090_2_0 n_20094_2_0 n_101160_2_0 
n_102720_2_0 n_101250_2_0 n_101260_2_0 n_101270_2_0 n_102800_2_0 n_102970_2_0; 
do i=1 to dim(s);
if n_20085_2_0 = 3 or   n_20085_2_0 = 4 or  n_20085_2_0 = 6  then s(i)=.;
end;
array t(*)  n_104000_3_0 n_104010_3_0 n_104020_3_0 n_104030_3_0 n_104040_3_0 n_104050_3_0 n_104060_3_0 n_104070_3_0 n_104080_3_0 n_104090_3_0 
n_104100_3_0 n_104110_3_0 n_104120_3_0 n_104130_3_0 n_104140_3_0 n_104150_3_0 n_104160_3_0 n_104170_3_0 n_104180_3_0 n_104190_3_0 n_104200_3_0 n_104210_3_0 
n_104220_3_0 n_104230_3_0 n_104240_3_0 n_104250_3_0 n_104260_3_0 n_104270_3_0 n_104280_3_0 n_104290_3_0 n_104300_3_0 n_104310_3_0 n_104320_3_0 n_104330_3_0 
n_104340_3_0 n_104350_3_0 n_104360_3_0 n_104370_3_0 n_104380_3_0 n_104400_3_0 n_104410_3_0 n_104420_3_0 n_104430_3_0 n_104440_3_0 n_104450_3_0 n_104460_3_0 
n_104470_3_0 n_104480_3_0 n_104490_3_0 n_104500_3_0 n_104510_3_0 n_104520_3_0 n_104530_3_0 n_104540_3_0 n_104550_3_0 n_104560_3_0 n_104570_3_0 n_104580_3_0 
n_104590_3_0 n_102410_3_0 n_102420_3_0 n_102430_3_0 n_102440_3_0 n_102450_3_0 n_103270_3_0 n_104000_3_0 n_104010_3_0 n_104110_3_0 n_104120_3_0 n_104280_3_0 
 n_100770_3_0 n_100800_3_0 n_100810_3_0 n_100820_3_0 n_100830_3_0 n_100840_3_0 n_100850_3_0 n_100860_3_0 n_102740_3_0 n_102780_3_0 n_102770_3_0 
n_102080_3_0 n_102090_3_0 n_20106_3_0  n_102810_3_0 n_102850_3_0 n_102870_3_0 n_100520_3_0 n_100920_3_0 n_100160_3_0 n_100170_3_0 n_100180_3_0 
n_100190_3_0 n_100200_3_0 n_100210_3_0 n_100220_3_0 n_100230_3_0 n_100370_3_0 n_100380_3_0 n_100490_3_0 n_100500_3_0 n_100540_3_0 n_100550_3_0 n_103000_3_0 
n_103010_3_0 n_103020_3_0 n_103030_3_0 n_103040_3_0 n_103070_3_0 n_103080_3_0 n_103090_3_0 n_30530_3_0 
n_103990_3_0 n_104400_3_0 n_102400_3_0 n_103310_3_0	n_20088_3_0	 n_103280_3_0 n_103260_3_0 n_103290_3_0 n_100760_3_0
n_100940_3_0 n_20091_3_0  n_100950_3_0 n_20092_3_0	n_101020_3_0 n_20093_3_0 n_101090_3_0 n_20094_3_0 n_101160_3_0 
n_102720_3_0 n_101250_3_0 n_101260_3_0 n_101270_3_0 n_102800_3_0 n_102970_3_0; 
do i=1 to dim(t);
if n_20085_3_0 = 3 or   n_20085_3_0 = 4 or  n_20085_3_0 = 6  then t(i)=.;
end;
array fh(*) n_104000_4_0 n_104010_4_0 n_104020_4_0 n_104030_4_0 n_104040_4_0 n_104050_4_0 n_104060_4_0 n_104070_4_0 n_104080_4_0 n_104090_4_0 
n_104100_4_0 n_104110_4_0 n_104120_4_0 n_104130_4_0 n_104140_4_0 n_104150_4_0 n_104160_4_0 n_104170_4_0 n_104180_4_0 n_104190_4_0 n_104200_4_0 n_104210_4_0 
n_104220_4_0 n_104230_4_0 n_104240_4_0 n_104250_4_0 n_104260_4_0 n_104270_4_0 n_104280_4_0 n_104290_4_0 n_104300_4_0 n_104310_4_0 n_104320_4_0 n_104330_4_0 
n_104340_4_0 n_104350_4_0 n_104360_4_0 n_104370_4_0 n_104380_4_0 n_104400_4_0 n_104410_4_0 n_104420_4_0 n_104430_4_0 n_104440_4_0 n_104450_4_0 n_104460_4_0 
n_104470_4_0 n_104480_4_0 n_104490_4_0 n_104500_4_0 n_104510_4_0 n_104520_4_0 n_104530_4_0 n_104540_4_0 n_104550_4_0 n_104560_4_0 n_104570_4_0 n_104580_4_0 
n_104590_4_0 n_102410_4_0 n_102420_4_0 n_102430_4_0 n_102440_4_0 n_102450_4_0 n_103270_4_0 n_104000_4_0 n_104010_4_0 n_104110_4_0 n_104120_4_0 n_104280_4_0 
 n_100770_4_0 n_100800_4_0 n_100810_4_0 n_100820_4_0 n_100830_4_0 n_100840_4_0 n_100850_4_0 n_100860_4_0 n_102740_4_0 n_102780_4_0 n_102770_4_0 
n_102080_4_0 n_102090_4_0 n_20106_4_0  n_102810_4_0 n_102850_4_0 n_102870_4_0 n_100520_4_0 n_100920_4_0 n_100160_4_0 n_100170_4_0 n_100180_4_0 
n_100190_4_0 n_100200_4_0 n_100210_4_0 n_100220_4_0 n_100230_4_0 n_100370_4_0 n_100380_4_0 n_100490_4_0 n_100500_4_0 n_100540_4_0 n_100550_4_0 n_103000_4_0 
n_103010_4_0 n_103020_4_0 n_103030_4_0 n_103040_4_0 n_103070_4_0 n_103080_4_0 n_103090_4_0 n_30530_4_0 n_103990_4_0 n_104400_4_0 n_102400_4_0 n_103310_4_0	n_20088_4_0	 n_103280_4_0 n_103260_4_0 n_103290_4_0 n_100760_4_0
n_100940_4_0 n_20091_4_0  n_100950_4_0 n_20092_4_0	n_101020_4_0 n_20093_4_0 n_101090_4_0 n_20094_4_0 n_101160_4_0 
n_102720_4_0 n_101250_4_0 n_101260_4_0 n_101270_4_0 n_102800_4_0 n_102970_4_0; 
do i=1 to dim(fh);
if n_20085_4_0 = 3 or   n_20085_4_0 = 4 or  n_20085_4_0 = 6  then fh(i)=.;
end;
run;

/*计算多轮调查的均数*/
data a.diet24opty1;
set diet24opty;
/*energy*/
energy_total=mean (of n_100002_0_0 n_100002_1_0 n_100002_2_0 n_100002_3_0 n_100002_4_0);
/*healthy*/
/*Vegetables: All vegetables except potatoes and legumes*/
mixveg=mean(of n_104060_0_0	n_104060_1_0 n_104060_2_0 n_104060_3_0 n_104060_4_0) ;
vegpiece=mean(of n_104070_0_0 n_104070_1_0 n_104070_2_0 n_104070_3_0 n_104070_4_0);
Coleslaw=mean(of n_104080_0_0 n_104080_1_0 n_104080_2_0 n_104080_3_0 n_104080_4_0);
Side=mean(of n_104090_0_0 n_104090_1_0 n_104090_2_0 n_104090_3_0 n_104090_4_0);
Avocado=mean(of n_104100_0_0 n_104100_1_0 n_104100_2_0 n_104100_3_0 n_104100_4_0);
Beetroot=mean(of n_104130_0_0 n_104130_1_0 n_104130_2_0 n_104130_3_0 n_104130_4_0);
Broccoli=mean(of n_104140_0_0 n_104140_1_0 n_104140_2_0 n_104140_3_0 n_104140_4_0);
Butternut=mean(of n_104150_0_0 n_104150_1_0 n_104150_2_0 n_104150_3_0 n_104150_4_0);
Cabbage=mean(of n_104160_0_0 n_104160_1_0 n_104160_2_0 n_104160_3_0 n_104160_4_0);
Carrot=mean(of n_104170_0_0 n_104170_1_0 n_104170_2_0 n_104170_3_0 n_104170_4_0);
Cauliflower=mean(of n_104180_0_0 n_104180_1_0 n_104180_2_0 n_104180_3_0 n_104180_4_0);
Celery=mean(of n_104190_0_0 n_104190_1_0 n_104190_2_0 n_104190_3_0 n_104190_4_0);
Courgette=mean(of n_104200_0_0 n_104200_1_0 n_104200_2_0 n_104200_3_0 n_104200_4_0);
Cucumber=mean(of n_104210_0_0 n_104210_1_0 n_104210_2_0 n_104210_3_0 n_104210_4_0);
Garlic=mean(of n_104220_0_0 n_104220_1_0 n_104220_2_0 n_104220_3_0 n_104220_4_0);
Leek=mean(of n_104230_0_0 n_104230_1_0 n_104230_2_0 n_104230_3_0 n_104230_4_0);
Lettuce=mean(of n_104240_0_0 n_104240_1_0 n_104240_2_0 n_104240_3_0 n_104240_4_0);
Mushroom=mean(of n_104250_0_0 n_104250_1_0 n_104250_2_0 n_104250_3_0 n_104250_4_0);
Onion=mean(of n_104260_0_0 n_104260_1_0 n_104260_2_0 n_104260_3_0 n_104260_4_0);
Parsnip=mean(of n_104270_0_0 n_104270_1_0 n_104270_2_0 n_104270_3_0 n_104270_4_0);
pepper=mean(of n_104290_0_0 n_104290_1_0 n_104290_2_0 n_104290_3_0 n_104290_4_0);
Spinach=mean(of n_104300_0_0 n_104300_1_0 n_104300_2_0 n_104300_3_0 n_104300_4_0);
Sprouts=mean(of n_104310_0_0 n_104310_1_0 n_104310_2_0 n_104310_3_0 n_104310_4_0);
Sweetcorn=mean(of n_104320_0_0 n_104320_1_0 n_104320_2_0 n_104320_3_0 n_104320_4_0);
Freshtomato=mean(of n_104340_0_0 n_104340_1_0 n_104340_2_0 n_104340_3_0 n_104340_4_0);
Tinnedtomato=mean(of n_104350_0_0 n_104350_1_0 n_104350_2_0 n_104350_3_0 n_104350_4_0);
Turnip=mean(of n_104360_0_0 n_104360_1_0 n_104360_2_0 n_104360_3_0 n_104360_4_0);
Watercress=mean(of n_104370_0_0 n_104370_1_0 n_104370_2_0 n_104370_3_0 n_104370_4_0);
Otherveg=mean(of n_104380_0_0 n_104380_1_0 n_104380_2_0 n_104380_3_0 n_104380_4_0);
/*Vegetables小计*/
Vegetables= sum (of mixveg vegpiece Coleslaw Side Avocado Beetroot Broccoli Butternut Cabbage Carrot Cauliflower Celery Courgette Cucumber Garlic Leek Lettuce 
Mushroom Onion Parsnip pepper Spinach Sprouts Sweetcorn Freshtomato Tinnedtomato Turnip Watercress Otherveg); 
if (n_103990_0_0^=. | n_103990_1_0^=. | n_103990_2_0^=.  | n_103990_3_0^=.  | n_103990_4_0^=.) & Vegetables=. then Vegetablesnew=0;else Vegetablesnew=Vegetables;
/*Fruit: All fruits and fruit juices*/
Stewedfruit=mean(of n_104410_0_0 n_104410_1_0 n_104410_2_0 n_104410_3_0 n_104410_4_0);
Prune=mean(of n_104420_0_0 n_104420_1_0 n_104420_2_0 n_104420_3_0 n_104420_4_0);
Dried=mean(of n_104430_0_0 n_104430_1_0 n_104430_2_0 n_104430_3_0 n_104430_4_0);
Mixedfruit=mean(of n_104440_0_0 n_104440_1_0 n_104440_2_0 n_104440_3_0 n_104440_4_0);
Apple=mean(of n_104450_0_0 n_104450_1_0 n_104450_2_0 n_104450_3_0 n_104450_4_0);
Banana=mean(of n_104460_0_0 n_104460_1_0 n_104460_2_0 n_104460_3_0 n_104460_4_0);
Berry=mean(of n_104470_0_0 n_104470_1_0 n_104470_2_0 n_104470_3_0 n_104470_4_0);
Cherry=mean(of n_104480_0_0 n_104480_1_0 n_104480_2_0 n_104480_3_0 n_104480_4_0);
Grapefruit=mean(of n_104490_0_0 n_104490_1_0 n_104490_2_0 n_104490_3_0 n_104490_4_0);
Grape=mean(of n_104500_0_0 n_104500_1_0 n_104500_2_0 n_104500_3_0 n_104500_4_0);
Mango=mean(of n_104510_0_0 n_104510_1_0 n_104510_2_0 n_104510_3_0 n_104510_4_0);
Melon=mean(of n_104520_0_0 n_104520_1_0 n_104520_2_0 n_104520_3_0 n_104520_4_0);
Orange=mean(of n_104530_0_0 n_104530_1_0 n_104530_2_0 n_104530_3_0 n_104530_4_0);
Satsuma=mean(of n_104540_0_0 n_104540_1_0 n_104540_2_0 n_104540_3_0 n_104540_4_0);
Peach=mean(of n_104550_0_0 n_104550_1_0 n_104550_2_0 n_104550_3_0 n_104550_4_0);
Pear=mean(of n_104560_0_0 n_104560_1_0 n_104560_2_0 n_104560_3_0 n_104560_4_0);
Pineapple=mean(of n_104570_0_0 n_104570_1_0 n_104570_2_0 n_104570_3_0 n_104570_4_0);
Plum=mean(of n_104580_0_0 n_104580_1_0 n_104580_2_0 n_104580_3_0 n_104580_4_0);
Otherfruit=mean(of n_104590_0_0 n_104590_1_0 n_104590_2_0 n_104590_3_0 n_104590_4_0);
Orangejuice=mean(of n_100190_0_0 n_100190_1_0 n_100190_2_0 n_100190_3_0 n_100190_4_0);
Grapefruit=mean(of n_100200_0_0 n_100200_1_0 n_100200_2_0 n_100200_3_0 n_100200_4_0);
Purefruitjuice=mean(of n_100210_0_0 n_100210_1_0 n_100210_2_0 n_100210_3_0 n_100210_4_0);
/*Fruits小计*/
Fruits=sum (of Stewedfruit Prune Dried Mixedfruit Apple Banana Berry Cherry Grapefruit Grape Mango Melon Orange Satsuma Peach Pear Pineapple Plum Otherfruit Orangejuice Grapefruit Purefruitjuice);
if (n_104400_0_0^=. | n_104400_1_0^=. | n_104400_2_0^=.  | n_104400_3_0^=.  | n_104400_4_0^=.) & Fruits=. then Fruitsnew=0; else Fruitsnew=Fruits;
/*nuts and legumes: Nuts and peanut butter, dried beans, peas, tofu*/
Saltpeanut=mean(of n_102410_0_0 n_102410_1_0 n_102410_2_0 n_102410_3_0 n_102410_4_0);
Unsaltpeanut=mean(of n_102420_0_0 n_102420_1_0 n_102420_2_0 n_102420_3_0 n_102420_4_0);
Saltnut=mean(of n_102430_0_0 n_102430_1_0 n_102430_2_0 n_102430_3_0 n_102430_4_0);
Unsaltnut=mean(of n_102440_0_0 n_102440_1_0 n_102440_2_0 n_102440_3_0 n_102440_4_0);
Seeds=mean(of n_102450_0_0 n_102450_1_0 n_102450_2_0 n_102450_3_0 n_102450_4_0);
Tofu=mean(of n_103270_0_0 n_103270_1_0 n_103270_2_0 n_103270_3_0 n_103270_4_0);
Bakedbean=mean(of n_104000_0_0 n_104000_1_0 n_104000_2_0 n_104000_3_0 n_104000_4_0);
Pulses=mean(of n_104010_0_0 n_104010_1_0 n_104010_2_0 n_104010_3_0 n_104010_4_0);
Broadbean=mean(of n_104110_0_0 n_104110_1_0 n_104110_2_0 n_104110_3_0 n_104110_4_0);
Greenbean=mean(of n_104120_0_0 n_104120_1_0 n_104120_2_0 n_104120_3_0 n_104120_4_0);
Pea=mean(of n_104280_0_0 n_104280_1_0 n_104280_2_0 n_104280_3_0 n_104280_4_0);
/*nuts小计*/
nuts=sum (of Saltpeanut Unsaltpeanut Saltnut Unsaltnut Seeds Tofu Bakedbean Pulses Broadbean Greenbean Pea );
if (n_102400_0_0^=. | n_102400_1_0^=. | n_102400_2_0^=. | n_102400_3_0^=. | n_102400_4_0^=. ) & nuts=. then nutsnew=0; else  nutsnew=nuts;
/* whole grains: Brown rice, dark breads, cooked cereal, whole grain cereal, ther grains, popcorn, wheat germ, bran*/
Porridge=mean(of n_100770_0_0 n_100770_1_0 n_100770_2_0 n_100770_3_0 n_100770_4_0);
Muesli=mean(of n_100800_0_0 n_100800_1_0 n_100800_2_0 n_100800_3_0 n_100800_4_0);
Oat=mean(of n_100810_0_0 n_100810_1_0 n_100810_2_0 n_100810_3_0 n_100810_4_0);
Sweetcereal=mean(of n_100820_0_0 n_100820_1_0 n_100820_2_0 n_100820_3_0 n_100820_4_0);
Plaincereal=mean(of n_100830_0_0 n_100830_1_0 n_100830_2_0 n_100830_3_0 n_100830_4_0);
Brancereal=mean(of n_100840_0_0 n_100840_1_0 n_100840_2_0 n_100840_3_0 n_100840_4_0);
Wholewheat=mean(of n_100850_0_0 n_100850_1_0 n_100850_2_0 n_100850_3_0 n_100850_4_0);
Othercereal=mean(of n_100860_0_0 n_100860_1_0 n_100860_2_0 n_100860_3_0 n_100860_4_0);
Brownrice=mean(of n_102740_0_0 n_102740_1_0 n_102740_2_0 n_102740_3_0 n_102740_4_0);
Othergrain=mean(of n_102780_0_0 n_102780_1_0 n_102780_2_0 n_102780_3_0 n_102780_4_0);
Couscous=mean(of n_102770_0_0 n_102770_1_0 n_102770_2_0 n_102770_3_0 n_102770_4_0);
/*Bread*/
/*wholemeal slice bread*/
if n_20091_0_0=3 then Slicedbread_0_0=n_100950_0_0; else Slicedbread_0_0=0;
if n_20091_0_0=. then Slicedbread_0_0=.;
if n_20091_1_0=3 then Slicedbread_1_0=n_100950_1_0; else Slicedbread_1_0=0;
if n_20091_1_0=. then Slicedbread_1_0=.;
if n_20091_2_0=3 then Slicedbread_2_0=n_100950_2_0; else Slicedbread_2_0=0;
if n_20091_2_0=. then Slicedbread_2_0=.;
if n_20091_3_0=3 then Slicedbread_3_0=n_100950_3_0; else Slicedbread_3_0=0;
if n_20091_3_0=. then Slicedbread_3_0=.;
if n_20091_4_0=3 then Slicedbread_4_0=n_100950_4_0; else Slicedbread_4_0=0;
if n_20091_4_0=. then Slicedbread_4_0=.;
Slicedbread=mean(of Slicedbread_0_0 Slicedbread_1_0 Slicedbread_2_0 Slicedbread_3_0 Slicedbread_4_0);
/*wholemeal Baguette*/
if n_20092_0_0=3 then Baguette_0_0=n_101020_0_0; else Baguette_0_0=0;
if n_20092_0_0=. then Baguette_0_0=.;
if n_20092_1_0=3 then Baguette_1_0=n_101020_1_0; else Baguette_1_0=0;
if n_20092_1_0=. then Baguette_1_0=.;
if n_20092_2_0=3 then Baguette_2_0=n_101020_2_0; else Baguette_2_0=0;
if n_20092_2_0=. then Baguette_2_0=.;
if n_20092_3_0=3 then Baguette_3_0=n_101020_3_0; else Baguette_3_0=0;
if n_20092_3_0=. then Baguette_3_0=.;
if n_20092_4_0=3 then Baguette_4_0=n_101020_4_0; else Baguette_4_0=0;
if n_20092_4_0=. then Baguette_4_0=.;
Baguette=mean(of Baguette_0_0 Baguette_1_0 Baguette_2_0 Baguette_3_0 Baguette_4_0);
/*wholemeal Bap*/
if n_20093_0_0=3 then Bap_0_0=n_101090_0_0; else Bap_0_0=0;
if n_20093_0_0=. then Bap_0_0=.;
if n_20093_1_0=3 then Bap_1_0=n_101090_1_0; else Bap_1_0=0;
if n_20093_1_0=. then Bap_1_0=.;
if n_20093_2_0=3 then Bap_2_0=n_101090_2_0; else Bap_2_0=0;
if n_20093_2_0=. then Bap_2_0=.;
if n_20093_3_0=3 then Bap_3_0=n_101090_3_0; else Bap_3_0=0;
if n_20093_3_0=. then Bap_3_0=.;
if n_20093_4_0=3 then Bap_4_0=n_101090_4_0; else Bap_4_0=0;
if n_20093_4_0=. then Bap_4_0=.;
Bap=mean(of Bap_0_0 Bap_1_0 Bap_2_0 Bap_3_0 Bap_4_0);
/*whole Breadroll*/
if n_20094_0_0=3 then Breadroll_0_0=n_101160_0_0; else Breadroll_0_0=0;
if n_20094_0_0=. then Breadroll_0_0=.;
if n_20094_1_0=3 then Breadroll_1_0=n_101160_1_0; else Breadroll_1_0=0;
if n_20094_1_0=. then Breadroll_1_0=.;
if n_20094_2_0=3 then Breadroll_2_0=n_101160_2_0; else Breadroll_2_0=0;
if n_20094_2_0=. then Breadroll_2_0=.;
if n_20094_3_0=3 then Breadroll_3_0=n_101160_3_0; else Breadroll_3_0=0;
if n_20094_3_0=. then Breadroll_3_0=.;
if n_20094_4_0=3 then Breadroll_4_0=n_101160_4_0; else Breadroll_4_0=0;
if n_20094_4_0=. then Breadroll_4_0=.;
Breadroll=mean(of Breadroll_0_0 Breadroll_1_0 Breadroll_2_0 Breadroll_3_0 Breadroll_4_0);
Wholemealpasta =mean(of n_102720_0_0 n_102720_1_0 n_102720_2_0 n_102720_3_0 n_102720_4_0);
Crispbread=mean(of n_101250_0_0 n_101250_1_0 n_101250_2_0 n_101250_3_0 n_101250_4_0);
Oatcakes=mean(of n_101260_0_0 n_101260_1_0 n_101260_2_0 n_101260_3_0 n_101260_4_0);
Otherbread=mean(of n_101270_0_0 n_101270_1_0 n_101270_2_0 n_101270_3_0 n_101270_4_0);
/*grains 小计*/
grains=sum (of  Porridge Muesli Oat Sweetcereal Plaincereal Brancereal Wholewheat Othercereal Brownrice Othergrain Couscous Slicedbread Baguette Bap Breadroll Wholemealpasta  Crispbread Oatcakes Otherbread);
if (n_100760_0_0^=. | n_100760_1_0^=. | n_100760_2_0^=. | n_100760_3_0^=. | n_100760_4_0^=. | n_100940_0_0^=. | n_100940_1_0^=. | n_100940_2_0^=. | n_100940_3_0^=. | n_100940_4_0^=.  ) & grains=. then grainsnew=0;
else grainsnew=grains;
/* low-fat dairy products : Skim milk, yogurt, cottage cheese*/
/*如果 yogurt类型是full fat 则处理为0，milk同样*/
if n_20106_0_0=210 then yogurt_0_0=n_102090_0_0;
if n_20106_1_0=210 then yogurt_1_0=n_102090_1_0;
if n_20106_2_0=210 then yogurt_2_0=n_102090_2_0;
if n_20106_3_0=210 then yogurt_3_0=n_102090_3_0;
if n_20106_4_0=210 then yogurt_4_0=n_102090_4_0;
if n_20106_0_0=211 then yogurt_0_0=0;
if n_20106_1_0=211 then yogurt_1_0=0;
if n_20106_2_0=211 then yogurt_2_0=0;
if n_20106_3_0=211 then yogurt_3_0=0;
if n_20106_4_0=211 then yogurt_4_0=0;
Yogurt=mean(of yogurt_0_0 yogurt_1_0 yogurt_2_0 yogurt_3_0 yogurt_4_0);
if n_100920_0_0=2102 or n_100920_0_0=2103  then milk_0_0=n_100520_0_0; else milk_0_0=0;
if n_100920_1_0=2102 or n_100920_1_0=2103  then milk_1_0=n_100520_1_0; else milk_1_0=0;
if n_100920_2_0=2102 or n_100920_2_0=2103  then milk_2_0=n_100520_2_0; else milk_2_0=0;
if n_100920_3_0=2102 or n_100920_3_0=2103  then milk_3_0=n_100520_3_0; else milk_3_0=0;
if n_100920_4_0=2102 or n_100920_4_0=2103  then milk_4_0=n_100520_4_0; else milk_4_0=0;
if  n_100920_0_0=. then milk_0_0=.;
if  n_100920_1_0=. then milk_1_0=.;
if  n_100920_2_0=. then milk_2_0=.;
if  n_100920_3_0=. then milk_3_0=.;
if  n_100920_4_0=. then milk_4_0=.;
Milk=mean(of milk_0_0 milk_1_0 milk_2_0 milk_3_0 milk_4_0);
hardcheese=mean(of n_102810_0_0 n_102810_1_0 n_102810_2_0 n_102810_3_0 n_102810_4_0);
cheesespread=mean(of n_102850_0_0 n_102850_1_0 n_102850_2_0 n_102850_3_0 n_102850_4_0);
Cottagecheese=mean(of n_102870_0_0 n_102870_1_0 n_102870_2_0 n_102870_3_0 n_102870_4_0);
cheese= sum (of hardcheese cheesespread Cottagecheese);
if (n_102800_0_0^=. | n_102800_1_0^=. |n_102800_2_0^=. | n_102800_3_0^=. | n_102800_4_0^=.) & cheese=. then cheesenew=0;
else cheesenew=cheese;
/*low-fat dairy 小计*/
lowfatdairy= sum (of Milk Yogurt cheesenew );
/*unhealthy components*/
/*sugar-sweetened drinks: Carbonated and noncarbonated sweetened beverages*/
/*废弃Lowdrink=mean(of n_100160_0_0 n_100160_1_0 n_100160_2_0 n_100160_3_0 n_100160_4_0);*/
Fizzydrink=mean(of n_100170_0_0 n_100170_1_0 n_100170_2_0 n_100170_3_0 n_100170_4_0);
Squash=mean(of n_100180_0_0 n_100180_1_0 n_100180_2_0 n_100180_3_0 n_100180_4_0);
Fruitsmootie=mean(of n_100220_0_0 n_100220_1_0 n_100220_2_0 n_100220_3_0 n_100220_4_0);
/*废弃Dairysmootie=mean(of n_100230_0_0 n_100230_1_0 n_100230_2_0 n_100230_3_0 n_100230_4_0);
sugarcoffee=mean(of n_100370_0_0 n_100370_1_0 n_100370_2_0 n_100370_3_0 n_100370_4_0);
artisugarcoffee=mean(of n_100380_0_0 n_100380_1_0 n_100380_2_0 n_100380_3_0 n_100380_4_0);
sugartea=mean(of n_100490_0_0 n_100490_1_0 n_100490_2_0 n_100490_3_0 n_100490_4_0);
artisugartea=mean(of n_100500_0_0 n_100500_1_0 n_100500_2_0 n_100500_3_0 n_100500_4_0);
Lowhotchocolate=mean(of n_100540_0_0 n_100540_1_0 n_100540_2_0 n_100540_3_0 n_100540_4_0);
Hotchocolate=mean(of n_100550_0_0 n_100550_1_0 n_100550_2_0 n_100550_3_0 n_100550_4_0);*/
/*sugar 小计*/
sugarsweetened=sum (of Fizzydrink Squash  Fruitsmootie );
/*Typical diet yesterday用做标识变量*/
if (n_100020_0_0^=. | n_100020_1_0^=. |n_100020_2_0^=. | n_100020_3_0^=. | n_100020_4_0^=.) & sugarsweetened=. then sugarsweetenednew=0; 
else sugarsweetenednew=sugarsweetened;
/*red and processed meats: Beef, pork, lamb, deli meats, organ meats, hot dogs, bacon*/
Sausage=mean(of n_103010_0_0 n_103010_1_0 n_103010_2_0 n_103010_3_0 n_103010_4_0);
Beef=mean(of n_103020_0_0 n_103020_1_0 n_103020_2_0 n_103020_3_0 n_103020_4_0);
Pork=mean(of n_103030_0_0 n_103030_1_0 n_103030_2_0 n_103030_3_0 n_103030_4_0);
Lamb=mean(of n_103040_0_0 n_103040_1_0 n_103040_2_0 n_103040_3_0 n_103040_4_0);
Bacon=mean(of n_103070_0_0 n_103070_1_0 n_103070_2_0 n_103070_3_0 n_103070_4_0);
Ham=mean(of n_103080_0_0 n_103080_1_0 n_103080_2_0 n_103080_3_0 n_103080_4_0);
Liver=mean(of n_103090_0_0 n_103090_1_0 n_103090_2_0 n_103090_3_0 n_103090_4_0);
Scotchegg=mean(of n_102970_0_0 n_102970_1_0 n_102970_2_0 n_102970_3_0 n_102970_4_0);
/*redmeat 小计*/
redmeat=sum (of Sausage Beef Pork Lamb Bacon Ham Liver Scotchegg);
if (n_103000_0_0^=. | n_103000_1_0^=. |n_103000_2_0^=. | n_103000_3_0^=. | n_103000_4_0^=.) & redmeat=. then redmeatnew=0;
else redmeatnew=redmeat;
/*estimated 24-hour sodium excretion: Sum of sodium content of all foods in FFQ*/
Sodium=n_30530_0_0;
run;
data a.dash;
set a.diet24opty1(keep =n_eid energy_total Vegetablesnew Fruitsnew nutsnew grainsnew lowfatdairy sugarsweetenednew redmeatnew Sodium);
run;
PROC EXPORT DATA= A.DASH 
            OUTFILE= "D:\1科研项目\0王巍巍课题\00博士阶段\2-导师\UKB\Depression comorbidity\STAT\dash.dta" 
            DBMS=STATA REPLACE;
RUN;

/*STATA中处理分位数,计算dashscore*/
/*use "D:\1科研项目\0王巍巍课题\00博士阶段\2-导师\UKB\Depression comorbidity\STAT\dash.dta"
 xtile quinVegetables =Vegetablesnew, nq(5)
 xtile quinFruits =Fruitsnew, nq(5)
 xtile quinnuts =nutsnew, nq(5)
 xtile quingrains =grainsnew, nq(5)
 xtile quinlowfatdairy =lowfatdairy, nq(5)
 xtile quinsugarsweetened =sugarsweetenednew, nq(5)
 xtile quinredmeat =redmeatnew, nq(5)
 xtile quinSodium =Sodium, nq(5)
 /*反向编码
 gen quinsugar=6-quinsugarsweetened if  quinsugarsweetened !=.
 gen quinmeat=6-quinredmeat if  quinredmeat !=.
 gen quinSodium_new=6-quinSodium if  quinSodium !=.
 egen dashscore=rowtotal ( quinVegetables  quinFruits quinnuts quingrains quinlowfatdairy quinsugarsweetened quinredmeat  quinSodium ) 
 replace dashscore=.  if  Vegetablesnew==. | Fruitsnew==. | nutsnew==. | grainsnew==. | lowfatdairy==. | sugarsweetenednew==. | redmeatnew==. | Sodium==./*191971*/
 /*均不缺失
 count if quinVegetables!=. & quinFruits!=. & quinnuts!=. & quingrains!=.   & quinlowfatdairy !=.  & quinsugarsweetened !=. & quinredmeat !=. &  quinSodium !=. 
 count if Vegetablesnew!=. & Fruitsnew!=. & nutsnew!=. & grainsnew!=.   & lowfatdairy !=.  & sugarsweetenednew!=. & redmeatnew !=. &  Sodium !=.
 centile dashscore, centile(23, 24, 25, 26 ,27, 28, 47, 48, 49, 50, 51, 52,71, 72, 73, 74, 75,76, 93, 94, 95)
 _pctile dashscore, p(24 25 49 50 74 75 94 95)  
 return list
 /*按位置划分
 sort dashscore
 generate newdashscore = group(4)
 gen dashpts=0 if  dashscore>0 & dashscore<17
 replace dashpts=25 if  dashscore<21 &  dashscore>=17
 replace dashpts=50 if  dashscore<26 &  dashscore>=21
 replace dashpts=80 if  dashscore<31 &  dashscore>=26
 replace dashpts=100 if  dashscore!=. &  dashscore>=31
 save "D:\1科研项目\0王巍巍课题\00博士阶段\2-导师\UKB\Depression comorbidity\STAT\dashscore.dta",replace*/


/*DASH之外的其他要素,数据来源于生活行为因素和biomarker*/
data a.DPRSTRAJv6;
set a.DPRSTRAJv5;
/*1 DASH diet score 变量名为dashpts*/
/*2 Physical activity score*/
/*S elf-reported minutes of moderate or vigorous physical activity per week. 
1 minute of vigorous physical activity is equivalent to 2 minutes of moderate physical activity*/
physical_mins=n_22038_0_0+n_22039_0_0*2;
if physical_mins >=150 then PA_pts=100;
else if physical_mins >=120 and physical_mins<150 then PA_pts=90;
else if physical_mins >=90 and physical_mins<120 then PA_pts=80;
else if physical_mins >=60 and physical_mins<90 then PA_pts=60;
else if physical_mins >=30 and physical_mins<60 then PA_pts=40;
else if physical_mins >=1 and physical_mins<30 then PA_pts=20;
else if physical_mins =0  then PA_pts=0;
/*3 Tobacco/nicotine exposure score*/
if  n_2897_0_0=-3 or n_2897_0_0=-1 then n_2897_0_0=.;
if  n_1249_0_0=-3 then  n_1249_0_0=.;
if  n_1259_0_0=-3 then  n_1259_0_0=.;
if n_2897_0_0>0 and age >=n_2897_0_0 then smkquit_y=age-n_2897_0_0;
if smoking=0 then smoke_cat=100;
else if smoking=1 and smkquit_y>=5 then smoke_cat=75;
else if n_1249_0_0=3 then smoke_cat=75;
else if smoking=1 and smkquit_y>=1 and smkquit_y <5 then smoke_cat=50;
else if n_1249_0_0=2 then smoke_cat=50;
else if smoking=1 and smkquit_y <1 then smoke_cat=25;
else if smoking=2 then smoke_cat=0;
if (n_1259_0_0=1 or n_1259_0_0=2) and smoke_cat>0  then smoke_pts=smoke_cat-20;
else  smoke_pts=smoke_cat;
/*4 Sleep health score*/
/*Self-reported average hours of sleep per night n_1160_0_0*/
if n_1160_0_0 =-1 or n_1160_0_0 =-3 then n_1160_0_0 =.;
if n_1160_0_0 >=7 and n_1160_0_0<9 then sleep_pts=100;
else if n_1160_0_0 >=9 and n_1160_0_0<10 then sleep_pts=90;
else if n_1160_0_0 >=6 and n_1160_0_0<7 then sleep_pts=70;
else if (n_1160_0_0 >=5 and n_1160_0_0<6) or n_1160_0_0 >=10 then sleep_pts=40;
else if n_1160_0_0 >=4 and n_1160_0_0<5 then sleep_pts=20;
else if n_1160_0_0 >=0 and n_1160_0_0<4 then sleep_pts=0;
/*5 Body mass index score*/
if bmi_i>0   and bmi_i< 25 then bmi_pts=100; 
else if bmi_i>=25 and bmi_i< 30 then bmi_pts=70; 
else if bmi_i>=30 and bmi_i< 35 then bmi_pts=30;
else if bmi_i>=35 and bmi_i< 40 then bmi_pts=15;
else if bmi_i>=40 then bmi_pts=0;
/*6 Blood lipid score: non-HDL cholesterol was calculated as total cholesterol minusHDLcholesterol*/
/*Units of measurement are mmol/L, LE8 unit 为mg/dl, 1mmol/L=38.67mg/dL cite PMID: 31017618  JAMA cardiol,to 
convert to millimoles per liter, multiply by 0.0259*/
nonhdl=(chol_i-hdl_chol_i)*38.67;
if chol_i < hdl_chol_i then nonhdl=0;
if nonhdl>=0 and nonhdl<130 then nonhdl_cat=100;
else if nonhdl>=130 and nonhdl<160 then nonhdl_cat=60;
else if nonhdl>=160 and nonhdl<190 then nonhdl_cat=40;
else if nonhdl>=190 and nonhdl<220 then nonhdl_cat=20;
else if nonhdl>=220 then nonhdl_cat=0;
if lipid_drug=1 and nonhdl_cat>0 then nonhdl_pts=nonhdl_cat-20; 
else nonhdl_pts=nonhdl_cat;
/*7 Glucose score*/
/*NGSP = 0.0915*IFCC + 2.15 cite Clin Chem 2019 PMID: 30518660 */
/*UKB中使用IFCC mmol/mol 作为unit*/
hba1c_n=hba1c_i*0.0915+ 2.15;
if diabetes_his=0 and hba1c_n>0  and hba1c_n<5.7 then hba1c_pts=100;
else if diabetes_his=0 and hba1c_n>=5.7 and hba1c_n<6.5 then hba1c_pts=60;
else if hba1c_n>0  and hba1c_n<7 then hba1c_pts=40;
else if  hba1c_n>=7 and hba1c_n<8 then hba1c_pts=30;
else if  hba1c_n>=8 and hba1c_n<9 then hba1c_pts=20;
else if  hba1c_n>=9 and hba1c_n<10 then hba1c_pts=10;
else if  hba1c_n>=10 then hba1c_pts=0;
/*8 Blood pressure score*/
if sbp_i>=160 or dbp_i>=100  then BP_cat=0; 
else if (sbp_i>=140 and sbp_i<160) or ( dbp_i>=90 and dbp_i<100) then BP_cat=25;
else if (sbp_i>=130 and sbp_i<140) or (dbp_i>=80 and dbp_i<90) then BP_cat=50;
else if (sbp_i>=120 and sbp_i<130) and ( dbp_i>=0 and dbp_i<80) then BP_cat=75;
else if sbp_i>0 and sbp_i<120 and dbp_i>0 and dbp_i<80 then BP_cat=100;
if hpt_drug=1 and BP_cat>0 then BP_pts=BP_cat-20;else BP_pts=BP_cat;
LE8score=mean (of dashpts PA_pts smoke_pts sleep_pts bmi_pts nonhdl_pts  hba1c_pts BP_pts);
if dashpts=. | PA_pts=. | smoke_pts=. | sleep_pts=. | bmi_pts=. |
nonhdl_pts=. |  hba1c_pts=. | BP_pts=. then LE8score=.;
if LE8score <50 and LE8score>=0 then CVH_cat=0;
if LE8score <80 and LE8score>=50 then CVH_cat=1;
if LE8score>=80 then CVH_cat=2;
run;
/*LE8*/
data A.LE8;
set a.DPRSTRAJv6 (keep=n_eid  LE8score CVH_cat dashpts PA_pts smoke_pts sleep_pts bmi_pts nonhdl_pts  hba1c_pts BP_pts );
run;
