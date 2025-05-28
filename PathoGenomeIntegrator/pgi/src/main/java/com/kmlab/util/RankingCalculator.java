package com.kmlab.util;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class RankingCalculator {
    private static final Logger logger = LogManager.getLogger(RankingCalculator.class);

    /**
     * 根据输入的整数列表计算每个整数的排名。
     * 列表中不同的整数将按照出现的顺序进行排名，相同的整数会分配相同的排名。
     * 排名从1开始。
     *
     * @param numbers 整数列表，将根据此列表计算排名。
     * @return 返回一个映射，其中包含列表中每个整数与其排名的对应关系。
     */
    public static Map<Integer, Integer> intMapRank(List<Integer> numbers) {
        logger.info("生成列表中每个整数与其排名的对应关系");
        Collections.sort(numbers);
        int rank = 0;
        Integer previousValue = null;
        Map<Integer, Integer> rankMap = new LinkedHashMap<>();
        for (Integer number : numbers) {
            if (previousValue == null || !number.equals(previousValue)) {
                rank++;
            }
            rankMap.put(number, rank);
            previousValue = number;
        }
        return rankMap;
    }
}
