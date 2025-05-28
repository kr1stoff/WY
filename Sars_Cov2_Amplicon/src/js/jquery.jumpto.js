/* ===========================================================
 * jquery-jumpto.js v1
 * ===========================================================
 * Copyright 2013 Pete Rojwongsuriya.
 * http://www.thepetedesign.com
 *
 * Create a smooth jump to sub navigational sidebar
 * with one js call
 *
 * https://github.com/peachananr/jumpto
 *
 * ========================================================== */

!(function ($) {
  var defaults = {
    firstLevel: "> h2",
    secondLevel: false,
    innerWrapper: ".jumpto-block",
    offset: 400,
    animate: 1000,
    navContainer: false,
    anchorTopPadding: 20,
    showTitle: "Jump To",
    closeButton: true,
  };


  $.fn.jumpto = function (options) {
    var settings = $.extend({}, defaults, options),
      el = $(this),
      html = "",
      block = $(settings.innerWrapper),
      selectors = "",
      title = "",
      close = "";

    el.addClass("jumpto-cotainer");

    getTitle = function(_this){ 
//      let regex = /<img\s+src=\"src\/image\/(no|yes|warn)\.png\"\s+class=\"check_icon\">/;
      let regex = /<img\s+class=\"check_icon\"[^>]*>/m;
      let match = _this.html().match(regex);
      match = match?match[0]:'';
      return  _this.text() + match;
    };

    redrawMenu = function () {
      $(selectors.slice(0, -2)).each(function (index) {
        if (isScrolledIntoView($(this))) {
          $(".jumpto-subnav a")
            .removeClass("active")
            .parent()
            .find(" a[href='#" + $(this).attr("id") + "']")
            .addClass("active");

          if (
            $("a[href='#" + $(this).attr("id") + "']")
              .parent()
              .parent()
              .hasClass("jumpto-second")
          ) {
            $("a[href='#" + $(this).attr("id") + "']")
              .closest(".jumpto-second")
              .show();
          }
          if (
            $("a[href='#" + $(this).attr("id") + "']")
              .parent()
              .parent()
              .hasClass("jumpto-first")
          ) {
            $("a[href='#" + $(this).attr("id") + "']")
              .closest(".jumpto-first")
              .find(".jumpto-second")
              .hide();
          }
          if (
            $("a[href='#" + $(this).attr("id") + "']")
              .parent()
              .find(".jumpto-second")
          ) {
            $("a[href='#" + $(this).attr("id") + "']")
              .parent()
              .find(".jumpto-second")
              .show();
          }
        }
      });
      if ($(document).scrollTop() >= settings.offset) {
        $(".jumpto-subnav").removeClass("bottom").addClass("fixed");
      } else {
        $(".jumpto-subnav").removeClass("bottom fixed");
      }
      if ($(document).scrollTop() >= el.outerHeight(true)) {
        $(".jumpto-subnav").addClass("bottom fixed");
      }
    };

    block.find(settings.firstLevel).each(function (index) {
      var b = $(this),
        i = index,
        inner_html = "";

      if (b.parent().find(settings.secondLevel).length > 0) {
        inner_html += "<ul class='jumpto-second'>";
        b.parent()
          .find(settings.secondLevel)
          .each(function (index) {
            var id = "jumpto_" + i + "_" + index;
            $(this).attr("id", id);
            link_to = "<a href='#" + id + "'>" + getTitle($(this))+ "</a>";
            inner_html += "<li>" + link_to + "</li>";
            selectors += "#" + id + ", ";
          });
        inner_html += "</ul>";
        var id = "jumpto_" + i;
        b.attr("id", id);
        link_to =
          "<a href='#" +
          id +
          "'>" +
           getTitle(b) + 
          "<img style='position: absolute;top: 8px;right: 0px;height:25px;width:25px' src='data:image/png;base64,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' /></a>";
        selectors += "#" + id + ", ";
        html += "<li>" + link_to + inner_html + "</li>";
      } else {
        var id = "jumpto_" + i;
        link_to = "<a href='#" + id + "'>" + getTitle(b) + "</a>";
        b.attr("id", id);
        selectors += "#" + id + ", ";
        html += "<li>" + link_to + "</li>";
      }
    });
    if (settings.showTitle != false) {
      var title = "<div class='jumpto-title'>" + settings.showTitle + "</div>";
    }

    if (settings.closeButton != false) {
      var close =
        "<div class='jumpto-close'><a href='#' id='jumpto-close'>Close</a></div>";
    }
    if (settings.navContainer == false) {
      $(this).append(
        "<nav class='jumpto-subnav'>" +
          title +
          "<ul class='jumpto-first'>" +
          html +
          "</ul>" +
          close +
          "</nav>"
      );
    } else {
      $(settings.navContainer)
        .addClass("jumpto-subnav")
        .html(title + "<ul class='jumpto-first'>" + html + "</ul>" + close);
    }

    $(".jumpto-subnav a[href*=#]:not([href=#])").click(function () {
      if ($(this).parent().find("ul").length > 0) {
        $(this).parent().find("ul").toggle(200);
      }

      if (
        location.pathname.replace(/^\//, "") ==
          this.pathname.replace(/^\//, "") ||
        location.hostname == this.hostname
      ) {
        var target = $(this.hash),interval = '',count = 0
        target = target.length
          ? target
          : $("[name=" + this.hash.slice(1) + "]");
        if (target.length) {
          $("html,body").animate(
            {
              scrollTop: target.offset().top - settings.anchorTopPadding,
            },
            settings.animate,
            "swing"
          );
          interval = setInterval(function() {
            //获取目标元素距离浏览器顶端的距离
            console.log(111111111);
            let topLength = target.offset().top - document.documentElement.scrollTop
            if((topLength > settings.anchorTopPadding + 5 || topLength < settings.anchorTopPadding - 5) && count < 5){
              $("html,body").animate(
                {
                  scrollTop: target.offset().top - settings.anchorTopPadding,
                },
                settings.animate,
                "swing"
              );
              count ++
            } else {
              clearInterval(interval);
            }
          },settings.animate)
          return false;
        }
      }
    });

    $(window).scroll(function () {
      redrawMenu();
    });

    $(".jumpto-subnav #jumpto-close").click(function () {
      var btn = $(this);
      btn
        .parent()
        .parent()
        .find("> .jumpto-first")
        .slideToggle("slow", function () {
          if ($(this).is(":visible")) {
            btn.html("Close");
          } else {
            btn.html("Open");
          }
        });
      return false;
    });

    setInterval(function () {
      var track = [];
      $(selectors.slice(0, -2)).each(function (index) {
        track.push(isScrolledIntoView($(this)));
      });
      if ($.inArray(true, track) == -1) {
        $(".jumpto-subnav a").removeClass("active");
        $(".jumpto-subnav .jumpto-second").hide();
      }
    }, 500);
    function isScrolledIntoView(elem) {
      //var docViewTop = $(window).scrollTop();
      //var docViewBottom = docViewTop + ($(window).height() /4);
  
      var elemTop = $(elem).offset().top;
      var elemBottom = elemTop + $(elem).parent().height();
      //return ((elemBottom <= docViewBottom) && (elemTop >= docViewTop));
      //return (elemTop >= docViewTop);
      //return ((elemTop <= $(window).scrollTop() + ($(window).height() / 10)) && (elemBottom >= $(window).scrollTop()));
      return (
        elemTop <= $(window).scrollTop() + (settings.anchorTopPadding + 2) &&
        elemBottom >= $(window).scrollTop() - (settings.anchorTopPadding - 2)
      );
    }
  };
})(window.jQuery);
