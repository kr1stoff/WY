$(function () {
  $(".fold1").click(function () {
    $(".jumpto-subnav").animate({ marginLeft: "-300px" });
    $(".toggleNav").animate({ left: "20px" });
    $("#report_body").animate({ paddingLeft: "60px" },function(){$.fn.dataTable.tables( { visible: true, api: true } ).columns.adjust();});
    $(this).addClass("close");
    $(".fold2").removeClass("close");
    

  });
  $(".fold2").click(function () {
    $(".jumpto-subnav").animate({ marginLeft: "0px" });
    $(".toggleNav").animate({ left: "260px" });
    $("#report_body").animate({ paddingLeft: "300px" },function(){$.fn.dataTable.tables( { visible: true, api: true } ).columns.adjust();});
    $(this).addClass("close");
    $(".fold1").removeClass("close");

  });
});
