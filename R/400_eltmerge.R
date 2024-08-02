###############################################################################
#' ELT CEP return period calculator
#'
#' @description
#' A utility that makes a vector of CEP return periods
#' that correspond to each of \code{nevents},
#' assuming the model covers \code{nyears},
#' using the midpoint plotting position rule to determine the probabilities
#' (but note that other rules are possible).
#'
#' @inheritParams man
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
make_cep_rps=function(nevents,nyears){

#
# 1) cep probabilities using the mid-point plotting position rule
#
	p=(c(1:nevents)-0.5)/nevents

#
# 2) cep frequencies from probabilities
#
	f=(nevents/nyears)*(1-p)

#
# 3) cep rps from frequencies
#
	rp=1/f

	return(rp)

}
###############################################################################
#' ELT CEP return period to index calculator
#'
#' #' @description
#' A simple utility that converts an RP to an index in an event set
#' given the total number of events, and the number of years,
#' using the midpoint plotting position rule to determine the probabilities
#' (but note that other rules are possible).
#'
#' @inheritParams man
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
convertrp2index=function(rp,nevents,nyears){

#
# 1) calculate the corresponding probability in the model losses
#
	p=nyears/(rp*nevents)
#
# 2) find the index for the corresponding model event
#
	index=max(1,min(round((1-p)*nevents+0.5),nevents))

	return(index)
}
###############################################################################
#' Adjusts a Catastrophe Model ELT using Historical Event Losses
#'
#' @description
#' \code{eltmerge} takes a cat model ELT (event loss table)
#' and a historical loss ELT and adjusts
#' the short RP (return period) losses in the model ELT to better match the historical losses.
#'
#' @inheritParams man
#'
#' @section Returns:
#'
#' \code{eltmerge} returns a list containing the following:
#'
#' \itemize{
#' \item \code{newmodelw1: }{adjusted model losses, using a weight of 1}
#' \item \code{newmodelw0: }{adjusted model losses, using a weight of 0}
#' \item \code{newmodelww: }{adjusted model losses,
#' using the user-defined weight}
#' \item \code{top: }{The index for the highest adjusted model loss}
#' }
#'
#' The adjusted model losses are in the same event order
#' as the original model losses.
#'
#' @section Details:
#'
#' \itemize{
#' \item Imagine you have a cat model ELT and some historical losses
#' \item The historical losses give you estimates of the short RP losses
#' \item The historical losses and the model don't agree very well, and you
#' tend to believe the historical losses more than the model losses
#' \item And so you'd like to adjust the short RP model losses towards the historical
#' losses
#' \item There are many ways that this problem could be addressed.
#' It's not rocket science, but it's also a bit tricky to sort out the
#' details.
#' \code{eltmerge} implements a relatively simple non-parametric method for
#' doing that. It is the simplest method I could come up with that
#' seemed to make sense and do what one might want.
#' \item It return three different adjusted versions of the model ELT
#' (\code{w1}, \code{w0} and \code{ww}).
#' \item The adjusted ELTs contains the same events as the input model ELT,
#' in the same order,
#' but with different losses on the short RP events.
#' \item The \code{w1} method works by taking each model event,
#' calculating the RP,
#' calculating the historical loss for that RP, and changing the loss
#' of the model event to be equal to that historical loss.
#' I call this method "RP loss matching".
#' \item The RPs of historical and model losses are calculated using the
#' mid-point plotting position method, linearly interpolated.
#' The mid-point method gives an RP for the
#' largest historical event of \code{2 x nyhist} years. Note that various methods
#' are used in the scientific literature for calculating plotting positions
#' and return periods,
#' in addition to the mid-point
#' method, and each would give different results. I think the mid-point
#' method is the best of the simple methods for calculating
#' plotting positions.
#' \item The \code{w0} method moves the results from the \code{w1}
#' method slightly
#' back towards to the original model losses. It does this by making a weighted
#' average of the results from the \code{w0} method and the original
#' model losses.
#' The weights on the \code{w0} results vary by return period, from a value
#' of 1 on the shortest return periods, to a value of 0 on the longest
#' historical loss return period (and linear as a function of RP in between).
#' \item The \code{ww} method gives results in between the \code{w1} method
#' and the \code{w0}, according to the user-defined weight.
#' \item The way the method works, and the difference between the three methods,
#' is perhaps best understood just by looking at the plots that are produced
#' by the example given below, which uses some realistic example data.
#' \item Please don't use this as a black box. There is nothing
#' scientifically "correct"
#' about the results: it is just an algorithm that tries to do something
#' sensible to solve this problem, based on various assumptions.
#' }
#'
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @seealso
#'
#' \itemize{
#' \item \code{eltmerge_plot: }{plots the results from \code{eltmerge},
#' in a 4 panel plot}
#' \item \code{eltmerge_plotone: }{plots the results from \code{eltmerge},
#' in a 1 panel plot}
#' \item \code{eltmerge_stats: }{calculates various statistics from the results
#' from \code{eltmerge}}
#' }
#' These routines are used in the example below.

#' @example man/examples/400_eltmerge_example.R
#'
#' @name eltmerge
#' @export
eltmerge=function(nyhist,histloss,nymodel,modelloss,weightw){
#
# 1 intro
#
	nehist=length(histloss)
	nemodel=length(modelloss)
#
# 2 sort input losses, and prepare to unsort later
#
# -hist
	sorthist=order(histloss)
	unsorthist=order(sorthist)
	shistloss=histloss[sorthist]
# -model
	sortmodel=order(modelloss)
	unsortmodel=order(sortmodel)
	smodelloss=modelloss[sortmodel]
#
# 3 make rps using a function
#
	rphist=make_cep_rps(nehist,nyhist)
	rpmodel=make_cep_rps(nemodel,nymodel)
#
# 4 find the model event that has the same RP as the largest hist event,
#		and index it 'top'
#
	top=convertrp2index(rphist[nehist],nemodel,nymodel)
#
# 5 find the largest historical rp
#
	rptop=rpmodel[top]
#
# 6 loop thru model events up to top
# -for each one, replace the modelled loss with the hist loss that
# 	has the same rp
# -and then replaced the modelled loss with that averaged with the
# 	modelled loss using a weight
	newmodelw1=smodelloss
	newmodelw0=smodelloss
	newmodelww=smodelloss
	for (i in 1:top){

# -calculate the model rp for model event i
		rp=rpmodel[i]

# -calculate the hist probability for that rp
		p1=1-nyhist/(nehist*rp)

# -calculate the nearest historical indices
		jm=max(1,min(nehist,round(nehist*p1-0)))
		jp=max(1,min(nehist,round(nehist*p1+1)))

# -calculate the linear interpolation factor for historical losses
# -dealing with the case where the first model rp is lower than the
# 	first historical rp
		fact=1
		if(rp<rphist[1]){
			fact=1-rp/rphist[jp]
		} else {
			drp=(rphist[jp]-rphist[jm])
			if(drp!=0){
				fact=1-(rp-rphist[jm])/drp
			}
		}

# -make the interpolated loss from the two end-point historical losses
		if(rp<rphist[1]){
			newmodelw1[i]=fact*0+(1-fact)*shistloss[jp]
		} else {
			newmodelw1[i]=fact*shistloss[jm]+(1-fact)*shistloss[jp]
		}

# -also make weighted versions in which the weights depend on the RP
# -weight is the weight that applies to the adjusted model
# -the weights always start at 1 at zero RP
#
# -this version reduces to zero weight at the last point
		weight1=(rptop-rp)/rptop
		newmodelw0[i]=weight1*newmodelw1[i]+(1-weight1)*smodelloss[i]
# -this version reduces to weightw at the last point
		weight2=1+rp*(weightw-1)/rptop
		newmodelww[i]=weight2*newmodelw1[i]+(1-weight2)*smodelloss[i]
	}
# sort all the new models
	newmodelw1=newmodelw1[unsortmodel]
	newmodelw0=newmodelw0[unsortmodel]
	newmodelww=newmodelww[unsortmodel]

	return(	list(
					newmodelw1=newmodelw1,
					newmodelw0=newmodelw0,
					newmodelww=newmodelww,
					top=top))
}
###############################################################################
#' One panel Plotting function for \code{eltmerge}
#'
#' #' @description
#' Plots CEPs for the output ELTs from \code{eltmerge}
#'
#' @inheritParams man
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/400_eltmerge_example.R
#'
#' @export
eltmerge_plotone=function(nyhist,shistloss,nymodel,smodelloss,
								 snewmodelw1,snewmodelw0,snewmodelww,top,rpmin,rpmax,main){

#
# 1 intro
#
	nehist=length(shistloss)
	nemodel=length(smodelloss)
#
# 2 make the rps
#
	rphist=make_cep_rps(nehist,nyhist)
	rpmodel=make_cep_rps(nemodel,nymodel)
#
# 3 calculate the ymax value based on the rpmax value
#
	irpmaxm=convertrp2index(rpmax,nemodel,nymodel)
	irpminm=convertrp2index(rpmin,nemodel,nymodel)
	irpmaxh=convertrp2index(rpmax,nehist,nyhist)
	irpminh=convertrp2index(rpmin,nehist,nyhist)

#	cat("irpminm=",irpminm,"\n")
	ymin=smodelloss[irpminm]
	ymax=smodelloss[irpmaxm]
#
# 4 plot the model losses
#
	plot(rpmodel[irpminm:nemodel],smodelloss[irpminm:nemodel],type="l",
		main=main,
		xlab="CEP RP",
		ylab="Event Loss",
		ylim=c(0,ymax),xlim=c(0,rpmax),col="black",lwd=2)
#
# 5 plot the historical losses
#
	lines(rphist[irpminh:nehist],shistloss[irpminh:nehist],col="red",lwd=2)
#
# 6 put an X on the model loss at the largest historical RP
#
	points(rpmodel[top],smodelloss[top],pch=4,lwd=2,col="black")
#
# 7 plot the adjusted modelled losses, based on weight of 1
#
	lines(rpmodel[irpminm:nemodel],snewmodelw1[irpminm:nemodel],col="blue")
#
# 8 plot the adjusted modelled losses, based on weight of 0
#
	lines(rpmodel[irpminm:nemodel],snewmodelw0[irpminm:nemodel],col="purple")
#
# 9 plot the adjusted modelled losses, based on user defined weight
#
	lines(rpmodel[irpminm:nemodel],snewmodelww[irpminm:nemodel],col="green")
#
# 10 legend
#
	legend(0.5*rpmax,0.55*ymax,
		c("model",
			"historical",
			"blend (weight of 1)",
			"blend (weight of 0)",
			"blend (weight)"),
		col=c("black","red","blue","purple","green"),lty=1)
}
###############################################################################
#' #' @description
#' Plots CEPs for the output ELTs from \code{eltmerge}, a four panel plot.
#'
#' @inheritParams man
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/400_eltmerge_example.R
#'
#' @export
eltmerge_plot=function(nyhist,histloss,nymodel,modelloss,
							newmodelw1,newmodelw0,newmodelww,top,rpmin){
#
# 1 intro
#
	nehist=length(histloss)
	shistloss=sort(histloss)
	smodelloss=sort(modelloss)
	snewmodelw1=sort(newmodelw1)
	snewmodelw0=sort(newmodelw0)
	snewmodelww=sort(newmodelww)
#
# 2 plot the four panels with different ranges
#
par(mfrow=c(2,2))
histlossmax=shistloss[nehist]
#
# (a)
#
rpmax=1000
#ymax=20*histlossmax
eltmerge_plotone(nyhist,shistloss,nymodel,smodelloss,
				snewmodelw1,snewmodelw0,snewmodelww,top,rpmin,rpmax,"(a) RP 0-1000")
#
# (b)
#
rpmax=100
eltmerge_plotone(nyhist,shistloss,nymodel,smodelloss,
				snewmodelw1,snewmodelw0,snewmodelww,top,rpmin,rpmax,"(b) RP 0-100")
#
# (c)
#
rpmax=50
eltmerge_plotone(nyhist,shistloss,nymodel,smodelloss,
				snewmodelw1,snewmodelw0,snewmodelww,top,rpmin,rpmax,"(c) RP 0-50 ")
#
# (d)
#
rpmax=5
eltmerge_plotone(nyhist,shistloss,nymodel,smodelloss,
				snewmodelw1,snewmodelw0,snewmodelww,top,rpmin,rpmax,"(d) RP 0-5")

}
###############################################################################
#' Statistics on the output from \code{eltmerge}
#'
#' #' @description
#' Calculates various statistics from the output ELTs from \code{eltmerge}
#'
#' @inheritParams man
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @export
eltmerge_stats=function(nyhist,shistloss,nymodel,smodelloss,
	snewmodelw1,snewmodelw0,snewmodelww){
#
# 1 intro
#
nehist=length(shistloss)
nemodel=length(smodelloss)
#
# 2 historical loss information
#
cat("\n")
cat("1: historical losses:\n")
cat(" number of events=",nehist,"\n")
cat(" number of years=",nyhist,"\n")
AAEhist=nehist/nyhist
AELhist=mean(shistloss)
AALhist=mean(shistloss)*AAEhist
cat(" AAE=",AAEhist,"\n")
cat(" AEL=",AELhist,"\n")
cat(" AAL=",AALhist,"\n")
#
# 3 model loss information
#
cat("\n")
cat("2: model losses:\n")
cat(" number of events=",nemodel,"\n")
cat(" number of years=",nymodel,"\n")
AAEmodel=nemodel/nymodel
AELmodel=mean(smodelloss)
AALmodel=mean(smodelloss)*AAEmodel
cat(" AAE=",AAEmodel,"\n")
cat(" AEL=",AELmodel,"\n")
cat(" AAL=",AALmodel,"\n")
#
# 4 ratios of model to historical
#
cat("\n")
cat("3: model to historical ratios:\n")
cat(" AAE=",AAEmodel/AAEhist,"\n")
cat(" AEL=",AELmodel/AELhist,"\n")
cat(" AAL=",AALmodel/AALhist,"\n")
#
# 5 AALs for the adjusted models, vs historical
#
cat("\n")
AALw1=mean(snewmodelw1)*AAEmodel
AALw0=mean(snewmodelw0)*AAEmodel
AALww=mean(snewmodelww)*AAEmodel
cat("4: AALs vs history:\n")
cat(" model =",AALmodel/AALhist,"\n")
cat(" new w1=",AALw1/AALhist,"\n")
cat(" new w0=",AALw0/AALhist,"\n")
cat(" new ss=",AALww/AALhist,"\n")
#
}


