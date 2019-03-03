# sequence lengths where original rules are below target:
below.target <- c(10:100)[bounds$pa.0.0 < 0.925]
n.below.target <- length(below.target)

# plots:
cola <- "red"
colb <- "orange"
colc <- "blue"
col.targ <- "brown"
# specificity:
par(mar=c(bottom=3, left=4, top=0, right=0)+.1)
uppspec <- ceiling(max(c(bounds$pa.0,bounds$pb.0,
                         bounds$pc.0))*100)/100
lowspec <- floor(min(c(bounds$pa.0,bounds$pb.0,
                       bounds$pc.0))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowspec,uppspec), xlab="", ylab="specificity")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
abline(h=.925, col=col.targ)
points(x=10:100, y=bounds$pa.0, type="l", col=cola)
points(x=10:100, y=bounds$pb.0, type="l", col=colb)
points(x=10:100, y=bounds$pc.0, type="l", col=colc)
points(x=below.target, y=rep(.925,n.below.target), pch=19)
# sensitivity p=0.4:
par(mar=c(bottom=2, left=4, top=0, right=0)+.1)
uppsens..4 <- ceiling(max(1-c(bounds$pa..4,bounds$pb..4,
                              bounds$pc..4))*100)/100
lowsens..4 <- floor(min(1-c(bounds$pa..4,bounds$pb..4,
                            bounds$pc..4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens..4,uppsens..4), xlab="", 
     ylab="sensitivity shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
points(x=10:100, y=1-bounds$pa..4, type="l", col=cola)
points(x=10:100, y=1-bounds$pb..4, type="l", col=colb)
points(x=10:100, y=1-bounds$pc..4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa..4,1-bounds$pb..4,1-bounds$pc..4)[below.target-9])
# sensitivity p=0.8:
uppsens..8 <- ceiling(max(1-c(bounds$pa..8,bounds$pb..8,
                              bounds$pc..8))*100)/100
lowsens..8 <- floor(min(1-c(bounds$pa..8,bounds$pb..8,
                            bounds$pc..8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens..8,uppsens..8), xlab="", 
     ylab="sensitivity shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=1-bounds$pa..8, type="l", col=cola)
points(x=10:100, y=1-bounds$pb..8, type="l", col=colb)
points(x=10:100, y=1-bounds$pc..8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa..8,1-bounds$pb..8,1-bounds$pc..8)[below.target-9])
# sensitivity p=1.2:
uppsens.1.2 <- ceiling(max(1-c(bounds$pa.1.2,bounds$pb.1.2,
                               bounds$pc.1.2))*100)/100
lowsens.1.2 <- floor(min(1-c(bounds$pa.1.2,bounds$pb.1.2,
                             bounds$pc.1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens.1.2,uppsens.1.2), xlab="", 
     ylab="sensitivity shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=1-bounds$pa.1.2, type="l", col=cola)
points(x=10:100, y=1-bounds$pb.1.2, type="l", col=colb)
points(x=10:100, y=1-bounds$pc.1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa.1.2,1-bounds$pb.1.2,1-bounds$pc.1.2)[below.target-9])
# sensitivity p=2:
uppsens.2 <- ceiling(max(1-c(bounds$pa.2,bounds$pb.2,
                             bounds$pc.2))*100)/100
lowsens.2 <- floor(min(1-c(bounds$pa.2,bounds$pb.2,
                           bounds$pc.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens.2,uppsens.2), xlab="", 
     ylab="sensitivity shift 2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=1-bounds$pa.2, type="l", col=cola)
points(x=10:100, y=1-bounds$pb.2, type="l", col=colb)
points(x=10:100, y=1-bounds$pc.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa.2,1-bounds$pb.2,1-bounds$pc.2)[below.target-9])
# positive likelihood ratios p=0.4:
upplrpos..4 <- ceiling(max(c(bounds$lrposa..4,bounds$lrposb..4,
                             bounds$lrposc..4))*100)/100
lowlrpos..4 <- floor(min(c(bounds$lrposa..4,bounds$lrposb..4,
                           bounds$lrposc..4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos..4,upplrpos..4), xlab="", 
     ylab="positive lhood ratio shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa..4, type="l", col=cola)
points(x=10:100, y=bounds$lrposb..4, type="l", col=colb)
points(x=10:100, y=bounds$lrposc..4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa..4,bounds$lrposb..4,bounds$lrposc..4)[below.target-9])
# positive likelihood ratios p=0.8:
upplrpos..8 <- ceiling(max(c(bounds$lrposa..8,bounds$lrposb..8,
                             bounds$lrposc..8))*100)/100
lowlrpos..8 <- floor(min(c(bounds$lrposa..8,bounds$lrposb..8,
                           bounds$lrposc..8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos..8,upplrpos..8), xlab="", 
     ylab="positive lhood ratio shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa..8, type="l", col=cola)
points(x=10:100, y=bounds$lrposb..8, type="l", col=colb)
points(x=10:100, y=bounds$lrposc..8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa..8,bounds$lrposb..8,bounds$lrposc..8)[below.target-9])
# positive likelihood ratios p=1.2:
upplrpos.1.2 <- ceiling(max(c(bounds$lrposa.1.2,bounds$lrposb.1.2,
                              bounds$lrposc.1.2))*100)/100
lowlrpos.1.2 <- floor(min(c(bounds$lrposa.1.2,bounds$lrposb.1.2,
                            bounds$lrposc.1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos.1.2,upplrpos.1.2), xlab="", 
     ylab="positive lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa.1.2, type="l", col=cola)
points(x=10:100, y=bounds$lrposb.1.2, type="l", col=colb)
points(x=10:100, y=bounds$lrposc.1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa.1.2,bounds$lrposb.1.2,bounds$lrposc.1.2)[below.target-9])
# positive likelihood ratios p=2:
upplrpos.2 <- ceiling(max(c(bounds$lrposa.2,bounds$lrposb.2,
                            bounds$lrposc.2))*100)/100
lowlrpos.2 <- floor(min(c(bounds$lrposa.2,bounds$lrposb.2,
                          bounds$lrposc.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos.2,upplrpos.2), xlab="", 
     ylab="positive lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa.2, type="l", col=cola)
points(x=10:100, y=bounds$lrposb.2, type="l", col=colb)
points(x=10:100, y=bounds$lrposc.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa.2,bounds$lrposb.2,bounds$lrposc.2)[below.target-9])

# negative likelihood ratios p=0.4:
upplrneg..4 <- ceiling(max(c(bounds$lrnega..4,bounds$lrnegb..4,
                             bounds$lrnegc..4))*100)/100
lowlrneg..4 <- floor(min(c(bounds$lrnega..4,bounds$lrnegb..4,
                           bounds$lrnegc..4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg..4,upplrneg..4), xlab="", 
     ylab="negative lhood ratio shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega..4, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb..4, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc..4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega..4,bounds$lrnegb..4,bounds$lrnegc..4)[below.target-9])
# negative likelihood ratios p=0.8:
upplrneg..8 <- ceiling(max(c(bounds$lrnega..8,bounds$lrnegb..8,
                             bounds$lrnegc..8))*100)/100
lowlrneg..8 <- floor(min(c(bounds$lrnega..8,bounds$lrnegb..8,
                           bounds$lrnegc..8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg..8,upplrneg..8), xlab="", 
     ylab="negative lhood ratio shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega..8, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb..8, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc..8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega..8,bounds$lrnegb..8,bounds$lrnegc..8)[below.target-9])
# negative likelihood ratios p=1.2:
upplrneg.1.2 <- ceiling(max(c(bounds$lrnega.1.2,bounds$lrnegb.1.2,
                              bounds$lrnegc.1.2))*100)/100
lowlrneg.1.2 <- floor(min(c(bounds$lrnega.1.2,bounds$lrnegb.1.2,
                            bounds$lrnegc.1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg.1.2,upplrneg.1.2), xlab="", 
     ylab="negative lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega.1.2, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb.1.2, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc.1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega.1.2,bounds$lrnegb.1.2,bounds$lrnegc.1.2)[below.target-9])
# negative likelihood ratios p=2:
upplrneg.2 <- ceiling(max(c(bounds$lrnega.2,bounds$lrnegb.2,
                            bounds$lrnegc.2))*100)/100
lowlrneg.2 <- floor(min(c(bounds$lrnega.2,bounds$lrnegb.2,
                          bounds$lrnegc.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg.2,upplrneg.2), xlab="", 
     ylab="negative lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega.2, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb.2, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega.2,bounds$lrnegb.2,bounds$lrnegc.2)[below.target-9])
