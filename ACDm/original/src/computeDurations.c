#include "header.h" 

//START---computeDurationsSubSec----------------------------//
void computeDurationsSubSec(int *y,
		int *M,
		int *d,
		int *h,
		int *m,  //5
		double *s,
		int *yDur,
		int *MDur,
		int *dDur,
		int *hDur,  //10
		int *mDur,
		double *sDur,
		int *vol,
		double *price,
		int *volDur,  //15
		double *priceDur,
		int *Ntrans,  //vector of the number of transactions for each duration
		double *dur,
		int *Ntime,  //number of all transaction
		int *Ndur,  //20
		double *open,
		double *close,
		int *durType,
		int *zeroDurHandeling,
		double *priceChange,  //25
		int *cumVol){

	int i, j = 0, lastMonth = 0, lastDay = 0, lastYear = 0, tempNtrans = 1,  tempVol = 0;
	double tempTurnover = 0, lastPrice = 0, secOfDay = 0, lastSecOfDay = 0, tempDur = 0;

	if((*durType == 1) && (*zeroDurHandeling == 1)){ // trade durations with 0s removed/aggregated
		for (i = 0; i < *Ntime-1; i++) {
			if (!(lastYear == *(y + i) && lastMonth == *(M + i) && lastDay == *(d + i))){ //checks if not same date as last transaction
				while((*(h + i) * 3600 + *(m + i) * 60 + *(s + i)) <= *open)	i++;  //"removes" observations before opening time
				lastYear = *(y + i); lastMonth = *(M + i); lastDay = *(d + i);
				lastSecOfDay = *open;
				tempVol = 0; tempTurnover = 0; tempNtrans = 1;
			}
			secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
			if(secOfDay <= *close){
				if(secOfDay == *(h+i+1)*3600 + *(m+i+1)*60 + *(s+i+1)){ //next duration is zero
					tempVol += *(vol+i);
					tempTurnover += *(vol+i) * (*(price+i));
					tempNtrans++;
					tempDur = secOfDay - lastSecOfDay;
				}
				else if(((*(h+i-1)*3600 + *(m+i-1)*60 + *(s+i-1))) == secOfDay){ //this duration is zero but not next
					tempVol += (*(vol+i));
					tempTurnover += (*(vol+i)) * (*(price+i));
					//adds a row to the new duration object:
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(volDur+j) = tempVol;
					*(priceDur+j) = (tempTurnover/tempVol);
					*(dur+j) = tempDur;
					*(Ntrans+j) = tempNtrans;

					lastSecOfDay = secOfDay;
					j++;
					tempVol = 0;
					tempTurnover = 0;
					tempNtrans = 1;
				}
				else{ //neither the current duration or the next is zero
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(volDur+j) = *(vol + i);
					*(priceDur+j) = *(price + i);
					*(dur+j) = secOfDay - lastSecOfDay;
					*(Ntrans+j) = 1;

					lastSecOfDay = secOfDay;
					j++;
					tempNtrans = 1;
				}
			}
		}

		//for the last observation:
		secOfDay = *(h + *Ntime - 1) * 3600 + *(m + *Ntime - 1) * 60 + *(s+ *Ntime - 1);
		if (secOfDay <= *close && lastYear == *(y + i) && lastMonth == *(M + i) //only writes it if it not after close and not new date
				&& lastDay == *(d + i)) {
			if(tempNtrans > 1){ //zero duration
				tempVol += *(vol+*Ntime - 1);
				tempTurnover += *(vol+*Ntime - 1)*(*(price+*Ntime - 1));
				//adds a row to the new duration object:
				*(yDur+j) = *(y + *Ntime - 1);
				*(MDur+j) = *(M + *Ntime - 1);
				*(dDur+j) = *(d + *Ntime - 1);
				*(hDur+j) = *(h + *Ntime - 1);
				*(mDur+j) = *(m + *Ntime - 1);
				*(sDur+j) = *(s + *Ntime - 1);
				*(volDur+j) = tempVol;
				*(priceDur+j) = tempTurnover/tempVol;
				*(dur+j) = secOfDay - (*(h + *Ntime - tempNtrans - 1)*3600 + *(m + *Ntime - tempNtrans - 1)*60 + *(s + *Ntime - tempNtrans - 1));
				*(Ntrans+j) = tempNtrans;
			}
			else { //not zero duration
				if(secOfDay <= *close) {
					*(yDur+j) = *(y + *Ntime - 1);
					*(MDur+j) = *(M + *Ntime - 1);
					*(dDur+j) = *(d + *Ntime - 1);
					*(hDur+j) = *(h + *Ntime - 1);
					*(mDur+j) = *(m + *Ntime - 1);
					*(sDur+j) = *(s + *Ntime - 1);
					*(volDur+j) = *(vol + *Ntime - 1);
					*(priceDur+j) = *(price + *Ntime - 1);
					*(dur+j) = secOfDay - (*(h + *Ntime - 2) * 3600 + *(m + *Ntime - 2) * 60 + *(s + *Ntime - 2));
					*(Ntrans+j) = 1;
				}
			}
		} else j--;
	}
	else if((*durType == 1) && (*zeroDurHandeling == 0)){ // trade durations with 0s kept
		for (i = 0; i < *Ntime; i++) {
			if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
				while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
				lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
				lastSecOfDay = *open;
				tempNtrans = 1;
			}
			secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
			if(secOfDay <= *close){
				//adds a row to the new duration object:
				*(yDur+j) = *(y + i);
				*(MDur+j) = *(M + i);
				*(dDur+j) = *(d + i);
				*(hDur+j) = *(h + i);
				*(mDur+j) = *(m + i);
				*(sDur+j) = *(s + i);
				*(volDur+j) = *(vol + i);
				*(priceDur+j) = *(price + i);
				*(dur+j) = secOfDay - lastSecOfDay;
				lastSecOfDay = secOfDay;
				j++;
			}
		}
	}
	else if(*durType == 2){ // price durations
		if(*zeroDurHandeling == 1){ //0s should  be aggregated
			for (i = 0; i < *Ntime - 1; i++) {
				if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
					while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
					if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
						lastPrice = *(price+i-1); //if there were transaction the first second of the trade day, the last price during that second will be used
						lastSecOfDay = *open;
					}
					else{
						lastPrice = *(price+i); //if there were no transactions the first second of the day, the first transaction of the day will be used
						lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
						i++;
					}
					lastYear = *(y + i); lastMonth = *(M + i); lastDay = *(d + i);
					tempVol = 0; tempTurnover = 0; tempNtrans = 0;
				}
				secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

				if(secOfDay <= *close){
					tempVol += *(vol+i);
					tempNtrans++;

					if(fabs(*(price+i) - lastPrice) >= *priceChange && secOfDay != lastSecOfDay){ //price difference is big enough and duration is not zero
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(volDur+j) = tempVol;
						*(priceDur+j) = *(price+i);
						*(dur+j) = secOfDay - lastSecOfDay;
						*(Ntrans+j) = tempNtrans;

						lastSecOfDay = secOfDay;
						lastPrice = *(price+i);
						j++;
						tempVol = 0;
						tempTurnover = 0;
						tempNtrans = 0;
					}
				}
			}
			j--;
		}
		if(*zeroDurHandeling == 0){ //0s should not be aggregated
			for (i = 0; i < *Ntime - 1; i++) {
				if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
					while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
					if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
						lastPrice = *(price+i-1); //if there were transaction the first second of the trade day, the last price during that second will be used
						lastSecOfDay = *open;
					}
					else{
						lastPrice = *(price+i); //if there were no transactions the first second of the day, the first transaction of the day will be used
						lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
						i++;
					}
					lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
					tempVol = 0; tempTurnover = 0; tempNtrans = 0;
				}
				secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

				if(secOfDay <= *close){
					tempVol += *(vol+i);
					tempNtrans++;

					if(fabs(*(price+i) - lastPrice) >= *priceChange){ //price difference is big enough
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(volDur+j) = tempVol;
						*(priceDur+j) = *(price+i);
						*(dur+j) = secOfDay - lastSecOfDay;
						*(Ntrans+j) = tempNtrans;

						lastSecOfDay = secOfDay;
						lastPrice = *(price+i);
						j++;
						tempVol = 0;
						tempTurnover = 0;
						tempNtrans = 0;
					}
				}
			}
			j--;
		}
	}
	else if(*durType == 3){ // volume durations
		if(*zeroDurHandeling == 1){ //0s should  be aggregated
			for (i = 0; i < *Ntime - 1; i++) {
				if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
					while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
					if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
						//if there were transaction the first second of the trade day, the last price during that second will be used
						lastSecOfDay = *open;
					}
					else{
						//if there were no transactions the first second of the day, the first transaction of the day will be used
						lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
						i++;
					}
					lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
					tempVol = 0; tempTurnover = 0; tempNtrans = 0;
				}
				secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

				if(secOfDay <= *close){
					tempVol += *(vol+i);
					tempNtrans++;

					if(tempVol >= *cumVol && secOfDay != lastSecOfDay){ //Cumulated volume is large enough and duration is not zero
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(volDur+j) = tempVol;
						*(priceDur+j) = *(price+i);
						*(dur+j) = secOfDay - lastSecOfDay;
						*(Ntrans+j) = tempNtrans;

						lastSecOfDay = secOfDay;
						j++;
						tempVol = 0;
						tempTurnover = 0;
						tempNtrans = 0;
					}
				}
			}
			j--;
		} else if(*zeroDurHandeling == 0){ //0s should  not be aggregated
			for (i = 0; i < *Ntime  - 1; i++) {
				if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){
					while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) < *open)	i++;
					if((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open){
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) == *open)	i++;
						//if there were transaction the first second of the trade day, the last price during that second will be used
						lastSecOfDay = *open;
					}
					else{
						//if there were no transactions the first second of the day, the first transaction of the day will be used
						lastSecOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
						i++;
					}
					lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
					tempVol = 0; tempTurnover = 0; tempNtrans = 0;
				}
				secOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);

				if(secOfDay <= *close){
					tempVol += *(vol+i);
					tempNtrans++;

					if(tempVol >= *cumVol){ //Cumulated volume is large enough
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(volDur+j) = tempVol;
						*(priceDur+j) = *(price+i);
						*(dur+j) = secOfDay - lastSecOfDay;
						*(Ntrans+j) = tempNtrans;

						lastSecOfDay = secOfDay;
						j++;
						tempVol = 0;
						tempTurnover = 0;
						tempNtrans = 0;
					}
				}
			}
			j--;
		}
	}

	*Ndur = j + 1;
}
//END---computeDurationsSubSec----------------------------//

//START---computeDurationsShort----------------------------//
void computeDurationsShort(int *y,
		int *M,
		int *d,
		int *h,
		int *m, //5
		double *s,
		int *yDur,
		int *MDur,
		int *dDur,
		int *hDur, //10
		int *mDur,
		double *sDur,
		double *dur,
		int *Ndur,
		int *Ntrans, //15
		int *Ntime,
		int *open,
		int *close,
		int *zeroDurHandeling){ //19

	int i, j = 0, lastMonth = 0, lastDay = 0, lastYear = 0, tempNtrans = 1;
	double sekOfDay = 0, lastSekOfDay = 0, tempDur = 0;


	if(*zeroDurHandeling == 1){ //0s removed/aggregated
		for (i = 0; i < *Ntime-1; i++) {
			if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
				while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
				lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
				lastSekOfDay = *open;
				tempNtrans = 1;
			}
			sekOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
			if(sekOfDay <= *close){
				if(sekOfDay == *(h+i+1)*3600 + *(m+i+1)*60 + *(s+i+1)){ //next duration is zero
					tempNtrans++;
					tempDur = sekOfDay - lastSekOfDay;
				}
				else if(((*(h+i-1)*3600 + *(m+i-1)*60 + *(s+i-1))) == sekOfDay){ //this duration is zero but not next
					//adds a row to the new duration object:
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(dur+j) = tempDur;
					*(Ntrans+j) = tempNtrans;

					lastSekOfDay = sekOfDay;
					j++;
					tempNtrans = 1;
				}
				else{ //neither the current duration nor the next is zero
					*(yDur+j) = *(y + i);
					*(MDur+j) = *(M + i);
					*(dDur+j) = *(d + i);
					*(hDur+j) = *(h + i);
					*(mDur+j) = *(m + i);
					*(sDur+j) = *(s + i);
					*(dur+j) = sekOfDay - lastSekOfDay;
					*(Ntrans+j) = 1;

					lastSekOfDay = sekOfDay;
					j++;
				}
			}
		}

		//for the last observation:
		sekOfDay = *(h + *Ntime - 1) * 3600 + *(m + *Ntime - 1) * 60 + *(s + *Ntime - 1);
		if (sekOfDay <= *close && lastYear == *(y + i) && lastMonth == *(M + i) //only writes it if it not after close and not new date
				&& lastDay == *(d + i)) {
			if(tempNtrans > 1){ //zero duration
				//adds a row to the new duration object:
				*(yDur+j) = *(y +  *Ntime - 1);
				*(MDur+j) = *(M +  *Ntime - 1);
				*(dDur+j) = *(d +  *Ntime - 1);
				*(hDur+j) = *(h +  *Ntime - 1);
				*(mDur+j) = *(m +  *Ntime - 1);
				*(sDur+j) = *(s +  *Ntime - 1);
				*(dur+j) = sekOfDay - (*(h + *Ntime - tempNtrans - 1) * 3600 + *(m + *Ntime - tempNtrans - 1) * 60 + *(s + *Ntime - tempNtrans - 1));
				*(Ntrans+j) = tempNtrans;
			} else { //not zero duration
				*(yDur+j) = *(y +  *Ntime - 1);
				*(MDur+j) = *(M +  *Ntime - 1);
				*(dDur+j) = *(d +  *Ntime - 1);
				*(hDur+j) = *(h +  *Ntime - 1);
				*(mDur+j) = *(m +  *Ntime - 1);
				*(sDur+j) = *(s +  *Ntime - 1);
				*(dur+j) = sekOfDay - (*(h + *Ntime - tempNtrans - 1) * 3600 + *(m + *Ntime - tempNtrans - 1) * 60 + *(s + *Ntime - tempNtrans - 1));
				*(Ntrans+j) = 1;
			}
		} else j--;
	} else if(*zeroDurHandeling == 0){ //zeroes not removed
		for (i = 0; i < *Ntime; i++) {
					if(!(lastYear == *(y+i) && lastMonth == *(M+i) && lastDay == *(d+i))){ //if new day
						while((*(h+i) * 3600 + *(m+i) * 60 + *(s+i)) <= *open)	i++; //removes transactions before opening
						lastYear = *(y+i); lastMonth = *(M+i); lastDay = *(d+i);
						lastSekOfDay = *open;
						tempNtrans = 1;
					}
					sekOfDay = *(h+i) * 3600 + *(m+i) * 60 + *(s+i);
					if(sekOfDay <= *close){
						//adds a row to the new duration object:
						*(yDur+j) = *(y + i);
						*(MDur+j) = *(M + i);
						*(dDur+j) = *(d + i);
						*(hDur+j) = *(h + i);
						*(mDur+j) = *(m + i);
						*(sDur+j) = *(s + i);
						*(dur+j) = sekOfDay - lastSekOfDay;
						lastSekOfDay = sekOfDay;
						j++;
					}
		}
	}

	*Ndur = j + 1;
}
//END---computeDurationsShort----------------------------//
