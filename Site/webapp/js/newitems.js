/******************************************************************************
 * functions to flag sidebar list items that have been added since client's
 * last visit. State is maintained in a cookie.
 ******************************************************************************/
var readListCookieName = 'sbi_read';
var oldHeadingPadBot;
var listItems = new Array();

/*
 *  return associative array where key = list IDs of read items
 */
function getReadFromCookie() {
  var readMap = {};
  var cookie = wdk.api.getCookie(readListCookieName);

  if (cookie == null) return readMap;
  
  jQuery(cookie.split(',')).each(function(i, val){
    readMap[val] = 1;
  });

  return readMap;
}

/*
 *  For each sidebar <li> item, background color those that are not in the 
 *  cookie.
 *  Expected minimal DOM branch:
 *     <a class="heading" href="#">
 *       <div class="menu_lefttop_drop">
 *         <ul id='?'>
 *           <li id='?'></li>
 *         </ul>
 *       </div>
 */
function flagUnreadListItems() {
  var readMap = getReadFromCookie();
  var totalUnreadCount = 0;
  var open = new Array();
  
  jQuery('a.heading').each(function(j){
    
    var sectUnreadCount = 0

    var section = jQuery(this).next('div.menu_lefttop_drop:first');
    var display = section.css("display");
    
    jQuery(section).
      children('ul').children('li[id]').each(function(k){
        
        listItems.push(this.id);
                  
        if ( ! readMap[this.id]) {
          this.style.backgroundColor='#ffffa0';
          this.style.margin='2px';
          //this.style.paddingLeft='1px';
          sectUnreadCount++;
          totalUnreadCount++;
        }

    });
    
    if (sectUnreadCount > 0) {
    
      var label = "<p class='unreadlabel'>";
      
      if (display == "none") label = label + "expand for ";
      label = label + sectUnreadCount + " new item" +
          ((listItems.length > 1) ? "s" : "") + "</p>"

      oldHeadingPadBot = jQuery(this).css('padding-bottom');
      jQuery(this).css({'padding-bottom' : '8px'});
      jQuery(this).append(label);
    }

    // div content is visible, consider them read without requiring user interaction.
    if (display != "none") open.push(this);
    
  });
  
  jQuery(open).each(function(){ putReadInCookie(this); });
  //console.log('totalUnreadCount ' + totalUnreadCount);
}

/*
 *  To be called when clicking <a class="heading"> to expand a specific
 *  subsection.
 *  Create a new cookie having the list from the original cookie
 *  plus all the <li> items in the specific <div class="menu_lefttop_drop">.
 *  See flagUnreadListItems() for expected DOM structure.
 *  Give the cookie to the client.
 */
function putReadInCookie(headernode) {
  var newCookieVal = new Array();
  var readMap = getReadFromCookie();
  jQuery(headernode).next('div.menu_lefttop_drop:first').
    children('ul').children('li[id]').each(function(k){
       readMap[this.id] = 1;
  });

  var key;
  for(key in readMap) {
      if (key == null || key == "") continue;
      if (jQuery.inArray(key, listItems) < 0) continue;
      newCookieVal.push(key);
  }
  
  jQuery(headernode).children('p:first').remove();
  jQuery(headernode).css({'padding-bottom' : oldHeadingPadBot})

  // expiration a year from today
  wdk.api.storeIntelligentCookie(readListCookieName, newCookieVal, 365, '/', wdk.api.secondLevelDomain());
}