setwd("C:/Users/jiehu/Desktop")
pacman::p_load(readxl, tidyverse, lubridate, reshape2, ggplot2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 深圳急救中心数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2011-2019年：摘机时间 -> 挂机时间 -> 派车时间 -> 出车时间 -> 到达现场时间 -> 上车时间 -> 送达时间
#	呼叫原因 | 患者症状 | 疾病类型 | 急救措施 | 地址类型 | 接车地址经度 | 接车地址纬度
#	到达现场体温	上车体温	到达医院体温	到达现场血压	上车血压	到达医院血压	到达现场脉搏	上车脉搏	到达医院脉搏	到达现场呼吸	上车呼吸	到达医院呼吸	到达现场意识	上车意识	到达医院意识
# 2020-2023年：开始受理时刻 -> 结束受理时刻 -> 收到指令时刻 -> 驶向现场时刻 -> 到达现场时刻 -> 病人上车时刻 -> 到达医院时刻
#	呼救原因 | ?呼救类型? | 疾病类型 | *事故类型* | 地址类型 | 接车地址经度 | 接车地址纬度
dir="D:/projects/01大学/02科研论文/120深圳/120数据/清洗后数据/"
for (year in 2011:2023) {
	dat <- read_excel(paste0(dir, year, '.xlsx'))
	print(paste(year, nrow(dat), ncol(dat)))
	if (year %in% 2011:2019) {dat <- dat %>% rename(开始时刻=摘机时间)}
	if (year %in% 2020:2023) {dat <- dat %>% rename(开始时刻=开始受理时刻)}
	dat$时刻 <- as_datetime(dat$开始时刻); dat$日期 <- as.Date(dat$时刻); dat$钟点 <- format(dat$时刻, "%H:%M:%S")
	eval(parse(text=paste0('dat', year, '<- dat')))
} 
grep("时间|时刻", names(dat2011), value=T); grep("体温", names(dat2019), value=T)
dat2019$体温 <- as.numeric(dat2019$到达医院体温); dat2019$体温[dat2019$体温>42] <- NA
str(dat2019$体温)
plot(dat2019$日期, dat2019$体温)
dat2019 <- dat2019 %>% mutate(
	时刻 = as_datetime(开始时刻), 
	日期 = as.Date(时刻), 
	钟点 <- format(时刻, "%H:%M:%S")
)
hist(dat$日期, breaks="days")
ggplot(dat, aes(x=日期)) + geom_density(); 
dat1 <- subset(dat, 事故类型 !="非重大事故"); table(dat1$事故类型)
ggplot(dat1, aes(x=日期, color=事故类型)) + geom_density()
ggplot(dat1, aes(日期, counts)) + geom_col(aes(fill = 事故类型)) + coord_flip() +  scale_x_date(date_labels = "%b %Y")
grps= unique(dat$事故类型); grps[grps != "非重大事故"]
pdf(paste0(f1_label,'.',f2_label,'.pdf'))
par(mfrow=c(6,2), mai=c(0.5,1,0.5,0.5))
for (year in 2012:2023) {
	eval(parse(text=paste0('dat <- dat', year)))
	
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 深圳疾控中心数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 0212以来本地病例
dat <- read_excel("D:/projects/01大学/01重大项目/022深圳防疫策略/01市CDC数据/温莹——阳性和密接/0212以来所有本地汇总表.xlsx") %>%
	select(全库序号, 管辖区, 街道, 省标准, 病例分类, 末次深圳市核酸阴性采样时间, 新增时间) %>% 
	rename(阴性时间=末次深圳市核酸阴性采样时间, 阳性时间=新增时间) %>% 
	mutate (
		阳性时间 = as.Date(as.numeric(阳性时间), origin = "1899-12-30"), 
		阳性时间天 = as.numeric(阳性时间 - as.Date("2022-02-12")), #difftime(day1, day2, unit="days")
		阴性时间 = as.Date(as.numeric(gsub('\\..*', '', 阴性时间)), origin = "1899-12-30"),
		风险分区 = recode(省标准,
			"封控区筛查" = "封控", "集中隔离" = "封控", "密接筛查" = "封控", "居家隔离" = "封控", 
			"管控区筛查" = "封控", "江西协查" = "封控", 
			"重点人员" = "社会面", 
			"社区筛查" = "社会面", "主动检测" = "社会面", "主动就诊" = "社会面"
		)
	) %>%
	filter ( 阳性时间 < as.Date("2022-04-01") )
hist(dat$阳性时间, breaks="days")
table(dat$风险分区)
dat1 <- dat %>%	group_by(阳性时间, 风险分区) %>% summarize(阳性感染者=n()) %>% ungroup 
ggplot(dat1, aes(x=阳性时间, y=阳性感染者, group=风险分区, colour=风险分区)) + 
	geom_line(aes(linetype=风险分区)) +
	geom_point(cex=2) + xlab("") + theme_bw() + 
	geom_vline(xintercept = as.numeric(as.Date(c("2022-03-14", "2022-03-20"))), col=c("red","green"), lwd=2,lty=1:2)
ggplot(dat1, aes(x=阳性时间天, y=as.numeric(阳性感染者), group=风险分区, colour=风险分区)) + 
  geom_line() +
  scale_colour_discrete(guide = 'none') +
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label=风险分区), method = list(dl.combine("first.points", "last.points")), cex = 0.8) 
## 境外输入
dat <- read_excel("D:/projects/01南科大/01重大课题/022 深圳防疫策略/01数据/温莹——阳性和密接/输入病例报告卡2022.xlsx")
dat <- subset(dat, select=c("入境前居住或旅行的国家或地区", "年龄", "备注")) 
dat1 <- subset(dat, grepl("宝安机场|宝安国际机场|航班", 备注))
dat1$indate <- gsub('.*[，|。|；]([0-9|\\.|于|年月|日|号|-]{2,12}).*宝安机场.*', '\\1', dat1$备注) ## "-"必须放到最后
dat1$dxdate <- gsub('.*[，|。|；](.*[0-9])([^0-9]*(疾控|海关|复核|重采样).*阳性).*', '\\1', dat1$备注) 
dat1$indate <- gsub('月|日|\\.', '-', dat1$indate) ## 不用"."，要不然被当初小数点处理
dat1$dxdate <- gsub('月|日|\\.', '-', dat1$dxdate)
write.table(dat1[,c("入境前居住或旅行的国家或地区", "年龄", "indate", "dxdate", "备注")], file="airline.txt", na="NA", sep="\t", row.names=F, quote=F, append=F)
dat <- read.table(file="airline.txt", header=T, sep="\t", as.is=T)[,1:4]
dat$indate <- as.Date(as.character(dat$indate), "%m-%d")
dat$dxdate <- as.Date(as.character(dat$dxdate), "%m-%d")
dat$days <- as.numeric(difftime(dat$dxdate, dat$indate, unit="days"))
dat1 <- dat %>% group_by(region) %>% filter(n() >1) %>% ungroup() # filter(n() >1)
hist(dat1$indate, breaks="weeks")
ggplot(dat1, aes(indate, days, col=region)) + geom_point()
