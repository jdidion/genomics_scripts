import xml.dom.minidom

class Element(xml.dom.minidom.Element):
    """
    Patch for the Element class to produce XML output acceptable to Google Maps
    http://www.ninemoreminutes.com/2009/12/google-maps-with-python-and-kml/
    """ 
    def writexml(self, writer, indent='', addindent='', newl=''):
        """
        indent = current indentation
        addindent = indentation to add to higher levels
        newl = newline string
        """
        writer.write(indent+'<' + self.tagName)

        attrs = self._get_attributes()
        a_names = attrs.keys()
        a_names.sort()

        for a_name in a_names:
            writer.write(' %s="' % a_name)
            xml.dom.minidom._write_data(writer, attrs[a_name].value)
            writer.write('"')
        if self.childNodes:
            newl2 = newl
            if len(self.childNodes) == 1 and \
                self.childNodes[0].nodeType == xml.dom.minidom.Node.TEXT_NODE:
                indent, addindent, newl = '', '', ''
            writer.write('>%s'%(newl))
            for node in self.childNodes:
                node.writexml(writer,indent+addindent,addindent,newl)
            writer.write("%s</%s>%s" % (indent,self.tagName,newl2))
        else:
            writer.write('/>%s'%(newl))

def xal_address(fields):
    city = element('Locality', children={'LocalityName': fields['city']})
    if 'street' in fields:
        street = element('Thoroughfare', Type='Street', 
            children={'AddressLine': fields['street']})
        city.appendChild(street)
    country = element('Country')
    country.appendChild(element('CountryNameCode', fields['country'], 
        Scheme='iso.3166-2'))
    if 'state' in fields:
        state = element('AdministrativeArea', children={
            'AdministrativeAreaName': fields['state']})
        state.appendChild(city)
        country.appendChild(state)
    else:
        country.appendChild(city)
    return country

def element(name, value=None, children=None, attrs={}, doc=xml.dom.minidom.Document(), ns=None, 
            cdata=False, **kwargs):
    e = doc.createElement(name) if ns is None else doc.createElementNS(ns, name)
    if value:
        if not isinstance(value, xml.dom.minidom.Node): 
            value = text(doc, value, cdata)
        e.appendChild(value)
    if children:
        if isinstance(children, dict):
            for k,v in children.items():
                if isinstance(v, list):
                    for val in v: e.appendChild(element(k, val, doc=doc))
                else: 
                    e.appendChild(element(k, v, doc=doc))
        else:
            for c in children: e.appendChild(c)
    for k,v in dict(attrs, **kwargs).items(): e.setAttribute(k, v)
    return e

def text(doc, value, cdata=False):
    if not isinstance(value, basestring):
        value = str(value)
    if cdata:
        return doc.createCDATASection(value) 
    else:
        return doc.createTextNode(value)