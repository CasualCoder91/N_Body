#include "ProgressBar.h"
#include <iostream>

ProgressBar::ProgressBar(float minValue, float maxValue){
    if (minValue < 0)
        this->shift = abs(minValue);
	this->minValue = minValue + shift;
	this->maxValue = maxValue + shift;
    this->active = (minValue < maxValue);
}

ProgressBar::ProgressBar(float minValue, float maxValue, bool active){
    if (minValue < 0)
        this->shift = abs(minValue);
    this->minValue = minValue + shift;
    this->maxValue = maxValue + shift;
    this->active = active;
}

void ProgressBar::Print(){
    if (this->active) {
        progress = currentValue / (maxValue - minValue);
        std::cout << ansi::foreground_red << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        if(currentValue!= maxValue)
            std::cout << "] " << int(progress * 100.0) << " %\r" << ansi::reset;
        else
            std::cout << "] " << int(progress * 100.0) << " %\n" << ansi::reset;
    }
}

void ProgressBar::Update(float currentValue){
    this->active = true;
    if (currentValue == -1)
        this->active = false;
    this->currentValue = currentValue+shift;
}

float ProgressBar::GetProgress(){
    return progress;
}

void ProgressBar::SetActive(bool value){
    this->active = value;
}

bool ProgressBar::GetActive(){
    return this->active;
}
