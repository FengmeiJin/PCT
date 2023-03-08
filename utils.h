//
// Created by Fengmei Jin on 12/7/2022.
//

#ifndef PCT_UTILS_H
#define PCT_UTILS_H

#include <sys/timeb.h>

// format from 0.01 to 0-01
static string float2str(float value) {
    string str = to_string(value);
    size_t found = str.find_last_of('.');
    str.replace(found, 1, "-");
    found = str.find_last_not_of('0');
    if (str.at(found) == '-') {
        return str.substr(0, found + 2);
    }
    else {
        return str.substr(0, found + 1);
    }
}

static long long millisecond() {
    time_t now = time(nullptr);    // current seconds based on current system
    char *dt = ctime(&now);    // convert now to string form (local date and time)

    timeb t{};
    ftime(&t);
    return t.time * 1000 + t.millitm;
}

static unsigned long timestr_timestamp(const string &datetime = "1970.01.01 00:00:00") {
    if (datetime.length() < 19) {
        cout << "invalid string - cant convert to timestamp\n";
    }

    struct tm tm{};
    tm.tm_year = atoi(datetime.substr(0, 4).c_str()) - 1900;    // years since 1900
    tm.tm_mon = atoi(datetime.substr(5, 2).c_str()) - 1;    //  months since January [0-11]
    tm.tm_mday = atoi(datetime.substr(8, 2).c_str());       //  day of the month [1-31]
    tm.tm_hour = atoi(datetime.substr(11, 2).c_str());
    tm.tm_min = atoi(datetime.substr(14, 2).c_str());
    tm.tm_sec = atoi(datetime.substr(17, 2).c_str());
//    char buff[80];
//    strftime(buff, 80, "%Y.%m.%d %H:%M:%S", &tm);
//    std::cout << "should be: " << std::string(buff) << "\n";

    return mktime(&tm); // return unit: second
}

static int computeDayGap(TIMETYPE currT, TIMETYPE prevT) {
    // current date/time based on current system
    time_t time_curr = currT;
    tm *tm_curr = localtime(&time_curr);
    int year_curr = 1900 + tm_curr->tm_year;
    int yday_curr = tm_curr->tm_yday;   /* days since January 1 [0-365] */

    time_t time_prev = prevT;
    tm *tm_prev = localtime(&time_prev);
    int year_prev = 1900 + tm_prev->tm_year;
    int yday_prev = tm_prev->tm_yday;
    if(year_curr == year_prev) {
        return yday_curr - yday_prev;
    }
    else {
        return yday_curr + 1 + (year_prev % 4 == 0 ? 366 : 365) - yday_prev;
    }
}

static void split(const string &str, vector<string> &tokens, const string &delim = ",") {
    tokens.clear();
//    auto start = str.find_first_not_of(delim, 0);
//    auto position = str.find_first_of(delim, start);
//    while (position != string::npos || start != string::npos) {
//        tokens.emplace_back(move(str.substr(start, position - start)));        // [start, position)
//        // next token
//        start = str.find_first_not_of(delim, position);
//        position = str.find_first_of(delim, start);
//    }

    size_t start = 0;
    size_t position = str.find_first_of(delim, start);
    while (start <= position && (start <= str.size() || position <= str.size())) {
        if (start == position) {     // imply an empty token (but only the delimiter)
            tokens.emplace_back("");
            start++;
        }
        else {
            tokens.emplace_back(move(str.substr(start, position - start)));        // [start, position)
            start = position + delim.size();
        }
        if (position > str.size()) {
            break;
        }
        position = str.find_first_of(delim, start);
    }
}

static bool sortByTime(std::pair<IDTYPE, TIMETYPE> &p1, std::pair<IDTYPE, TIMETYPE> &p2) {
    return p1.second < p2.second;   // ascending order
}

static bool sortByTime2(std::tuple<TIMETYPE, float, float> &p1, std::tuple<TIMETYPE, float, float> &p2) {
    return std::get<0>(p1) < std::get<0>(p2);
}

static bool sortByCnt(pair<string, int> &p1, pair<string, int> &p2) {
    return p1.second > p2.second;   // descending
}

#endif //PCT_UTILS_H
